# =============================================================================
# Analysis 14 — Cross-fluid protein-protein correlation networks + Louvain
# community detection (inputs for Figure 4 d/e and Extended Data 1).
#
# For each timepoint, computes all pairwise Spearman correlations among the
# plasma, saliva and urine proteins, keeps edges passing an adjusted-p-value and
# bootstrap-stability threshold, and detects communities (modules) with
# igraph::cluster_louvain(). set.seed(1234) fixes the bootstrap; set.seed(123)
# fixes the network layout. Computationally heavy (all pairwise correlations
# across three fluids x timepoints, with bootstrap resampling).
#
# Input:  data/res.{saliva,urine,olink}.linear.rda, data/prot.label.{saliva,urine}.rda,
#         data/olink.prot.label.rda, data/exerome.dat.{saliva,urine}.rda,
#         data/olink.exerome.dat.rda        (from analyses 02 / 03 / 04)
# Output (data/network/):
#   network_crossfluid_padj_ORAaware_<tp>.graphml  (per-timepoint network: nodes/edges/modules/hubs)
#   ORA_sig_crossfluid_padj_<tp>_Module_<m>.csv    (which modules are ORA-significant)
#   network_metrics.csv                            (edge counts per timepoint)
# =============================================================================
# ============================== SETUP =======================================
suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(purrr)
    library(corrplot)
    library(igraph)
    library(ggplot2)
    library(doParallel)
    library(foreach)
    library(future)
    library(furrr)
})

set.seed(1234)  # reproducibility for bootstrap & layouts

# --------------------- Palettes & timepoint order ---------------------------
fluid_color_palette <- c(
    "plasma" = "#4B0082",
    "saliva" = "#FFA07A",
    "urine"  = "#FFD700"
)

time_color_palette <- c(
    "-1" = "#E69F00",
    "0"  = "#56B4E9",
    "0.5"= "#009E73",
    "1"  = "#F0E442",
    "3"  = "#0072B2",
    "24" = "#CC79A7"
)

# --------- Project root & output dirs (all within this repo) -----------------
PROJECT_ROOT <- tryCatch(here::here(), error = function(e) normalizePath(".", mustWork = FALSE))

OUT_DIR   <- PROJECT_ROOT
FILES_DIR <- file.path(PROJECT_ROOT, "data", "network")
FIG_DIR   <- file.path(PROJECT_ROOT, "figures", "crossfluid_networks")

dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(FILES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,   showWarnings = FALSE, recursive = TRUE)

# tiny helpers so we never hardcode strings again
p_files <- function(...) file.path(FILES_DIR, ...)
p_fig   <- function(...) file.path(FIG_DIR, ...)

# ------------------------ Analysis parameters -------------------------------
N_MIN         <- 16        # require >15 overlapping subjects per pair
R_MIN         <- 0.70      # absolute Spearman rho threshold for edges
ALPHA_ADJ     <- 0.05      # adjusted p-value threshold for edges
PADJ_METHOD   <- "fdr"
B_BOOT        <- 200       # bootstrap resamples per kept edge
STABILITY_MIN <- 0.70      # keep edges stable in >=70% of resamples

# Optional QA: re-check a random subset of correlations with base::cor() vs cor.test()
CHECK_COR     <- TRUE
CHECK_COR_N   <- 50        # max pairs to re-check per timepoint

# NEW: Parallel switches (set to FALSE to run sequentially)
PARALLEL_COR  <- TRUE      # parallelize pairwise cor.test loop
PARALLEL_BOOT <- TRUE      # parallelize bootstrap stability
N_WORKERS     <- max(1L, parallel::detectCores(logical = TRUE) - 1L)

# CROSS-FLUID switch: keep only cross-fluid edges if TRUE
CROSS_FLUID_ONLY <- TRUE

# ===== PDF/label settings to avoid encoding warnings =====
FLUID_SEP <- "-"            # use ASCII hyphen instead of en-dash in facet labels
USE_CAIRO <- FALSE          # set TRUE if you want Unicode via Cairo PDF (needs Cairo)

# Prevent BLAS over-subscription inside workers (harmless if no MKL/OpenBLAS)
Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
)

# Start parallel backends
.cl <- NULL
if (PARALLEL_COR) {
    .cl <- parallel::makeCluster(N_WORKERS)
    doParallel::registerDoParallel(.cl)
}
if (PARALLEL_BOOT) {
    future::plan(multisession, workers = N_WORKERS)   # use multicore on Linux/Mac if you prefer
}

# Which fluids are measured at each timepoint (hard guard)
allowed_fluids_by_tp <- list(
    "-1" = c("plasma","saliva","urine"),
    "0"  = c("plasma","saliva","urine"),
    "0.5"= c("plasma","saliva"),
    "1"  = c("plasma","saliva"),
    "3"  = c("plasma","saliva"),
    "24" = c("plasma","saliva","urine")
)

# ----------------------------- Helpers --------------------------------------
# Parse fluid from "<Assay>_<fluid>" safely (also tolerates ".1" suffixes)
get_fluid_from_feature <- function(s) {
    f <- sub(".*_", "", s)
    f <- tolower(trimws(f))
    f <- sub("\\..*$", "", f)
    f
}

# Linear rescale to 0..1 with a safe fallback (used for plot widths/sizes)
rescale01 <- function(x) {
    xr <- range(x, na.rm = TRUE)
    if (!is.finite(xr[1]) || !is.finite(xr[2]) || diff(xr) == 0) return(rep(0.5, length(x)))
    (x - xr[1]) / diff(xr)
}

# Normalize time labels like "0h", "24 h" -> "0", "24"
normalize_time <- function(x) {
    x <- trimws(as.character(x))
    x <- gsub("\\s*[hH]$", "", x)
    x
}

# Bootstrap stability for one pair (down-weights correlations driven by few subjects)
boot_stability_pair <- function(x, y, B = B_BOOT, r_min = R_MIN, n_min = N_MIN) {
    ok  <- is.finite(x) & is.finite(y)
    idx <- which(ok)
    if (length(idx) < n_min) return(0)
    x0 <- x[idx]; y0 <- y[idx]
    if (length(unique(x0)) < 2 || length(unique(y0)) < 2) return(0)
    rho0 <- suppressWarnings(cor(x0, y0, method = "spearman"))
    if (!is.finite(rho0)) return(0)
    s0 <- sign(rho0)
    pass <- 0L; n_ok <- length(idx)
    for (b in seq_len(B)) {
        idxb <- sample.int(n_ok, size = n_ok, replace = TRUE)
        xb <- x0[idxb]; yb <- y0[idxb]
        if (length(unique(xb)) < 2 || length(unique(yb)) < 2) next
        rhob <- suppressWarnings(cor(xb, yb, method = "spearman"))
        if (is.finite(rhob) && abs(rhob) >= r_min && sign(rhob) == s0) pass <- pass + 1L
    }
    pass / B
}

# ============================== LOAD DATA ===================================
load("data/res.urine.linear.rda")
load("data/prot.label.urine.rda")
load("data/exerome.dat.urine.rda")

load("data/res.saliva.linear.rda")
load("data/prot.label.saliva.rda")
load("data/exerome.dat.saliva.rda")

load("data/res.olink.linear.rda")
load("data/olink.exerome.dat.rda")
load("data/olink.prot.label.rda")

# ============================== TIDY HELPERS =================================
find_col <- function(df, candidates) {
    nms <- names(df)
    match_idx <- match(tolower(candidates), tolower(nms))
    col <- nms[na.omit(match_idx)][1]
    if (length(col) == 0) NA_character_ else col
}

add_id_and_drop <- function(df) {
    subj_col <- find_col(df, c("label","subject","id","sample","participant"))
    tf_col   <- find_col(df, c("t.factor","t_factor","time_label","timepoint","time","t"))
    if (is.na(subj_col) || is.na(tf_col)) {
        stop("Could not find subject/label and/or t.factor (time) columns in exerome.dat.")
    }
    df %>%
        mutate(
            .subject_clean = str_replace_all(.data[[subj_col]] %>% as.character(), "EX", ""),
            .t_char        = normalize_time(.data[[tf_col]]),
            ID             = paste0(.subject_clean, "_", .t_char)
        ) %>%
        select(-all_of(c(subj_col, tf_col, ".subject_clean", ".t_char")))
}

# Prepare per-fluid LONG dataframe of significant assays
prep_exeromedat <- function(exeromedat, prot_label, res_linear, fluid = "fluid") {
    if (anyDuplicated(names(exeromedat))) {
        warning("Duplicate column names in exeromedat; keeping first occurrence only.")
        exeromedat <- exeromedat[, !duplicated(names(exeromedat)), drop = FALSE]
    }
    if (anyDuplicated(names(prot_label))) {
        warning("Duplicate column names in prot_label; keeping first occurrence only.")
        prot_label <- prot_label[, !duplicated(names(prot_label)), drop = FALSE]
    }
    if (anyDuplicated(names(res_linear))) {
        res_linear <- res_linear[, !duplicated(names(res_linear)), drop = FALSE]
    }
    assay_col <- find_col(prot_label, c("Assay","assay"))
    if (is.na(assay_col)) stop("prot.label is missing an 'Assay' column.")
    if (assay_col != "Assay") prot_label <- prot_label %>% rename(Assay = all_of(assay_col))
    else                      prot_label <- prot_label %>% mutate(Assay = as.character(Assay))
    if (anyDuplicated(prot_label$Assay)) {
        warning("Duplicate Assay entries in prot_label; keeping first occurrence only.")
        prot_label <- prot_label %>% distinct(Assay, .keep_all = TRUE)
    }
    # Optional OlinkID -> Assay mapping on exeromedat
    olinkid_col <- find_col(prot_label, c("OlinkID","Olink.Id","Olink_ID","olinkid"))
    if (!is.na(olinkid_col)) {
        map_tbl <- prot_label %>%
            select(Assay, !!olinkid_col) %>%
            filter(!is.na(.data[[olinkid_col]]), .data[[olinkid_col]] != "", Assay != "") %>%
            distinct(.data[[olinkid_col]], .keep_all = TRUE)
        map_vec <- setNames(map_tbl$Assay, map_tbl[[olinkid_col]])
        overlap <- intersect(names(exeromedat), names(map_vec))
        if (length(overlap)) {
            names(exeromedat)[match(overlap, names(exeromedat))] <- unname(map_vec[overlap])
            exeromedat <- exeromedat[, !duplicated(names(exeromedat)), drop = FALSE]
        }
    }
    exeromedat2 <- add_id_and_drop(exeromedat)
    exeromedat2 <- exeromedat2[, !duplicated(names(exeromedat2)), drop = FALSE]
    keep_assays <- intersect(names(exeromedat2), prot_label$Assay)
    exeromedat2 <- exeromedat2 %>% select(ID, all_of(keep_assays))
    res_assay_col <- find_col(res_linear, c("Assay","assay","feature","gene","protein"))
    fdr_col       <- find_col(res_linear, c("fdr.aov","fdr","padj","adj.p","qvalue"))
    if (is.na(res_assay_col) || is.na(fdr_col)) {
        stop("res.linear is missing 'Assay' and/or FDR-like column.")
    }
    sig_assays <- res_linear %>%
        mutate(Assay = .data[[res_assay_col]] %>% as.character(),
               fdr   = as.numeric(.data[[fdr_col]])) %>%
        filter(!is.na(fdr), fdr < 0.05, Assay != "") %>%
        distinct(Assay, .keep_all = TRUE) %>%
        pull(Assay)
    keep_sig <- intersect(keep_assays, sig_assays)
    if (length(keep_sig) == 0) {
        warning(sprintf("No significant assays (fdr<0.05) for '%s'. Returning empty set.", fluid))
    }
    exerome_sig <- exeromedat2 %>% select(ID, all_of(keep_sig))
    long <- exerome_sig %>%
        pivot_longer(-ID, names_to = "Assay", values_to = "value") %>%
        distinct(ID, Assay, .keep_all = TRUE) %>%
        mutate(fluid = fluid)
    long
}

# ========================= BUILD LONG + WIDE ================================
urine_long <- prep_exeromedat(exerome.dat.urine,  prot.label.urine,  res.urine.linear,  "urine")
saliva_long<- prep_exeromedat(exerome.dat.saliva, prot.label.saliva, res.saliva.linear, "saliva")
olink_long <- prep_exeromedat(olink.exerome.dat,  olink.prot.label,  res.olink.linear,  "plasma")

combined_exeromedat_long <- bind_rows(urine_long, saliva_long, olink_long)
# combined_exeromedat_long: ID, Assay, value, fluid

# =========================== CORRELATIONS ===================================
# Tidy + numeric coercion (CRITICAL)
prep_long <- function(df_long) {
    df_long %>%
        mutate(
            subject = sub("_(?=[^_]+$)", "", ID, perl = TRUE),  # everything before LAST '_'
            time    = sub("^.*_", "", ID),
            feature = paste0(Assay, "_", fluid),
            value   = suppressWarnings(as.numeric(trimws(as.character(value))))
        ) %>%
        filter(is.finite(value)) %>%
        select(subject, time, feature, value, fluid)
}

# One timepoint -> wide matrix (rows=subjects, cols=features)
tp_to_matrix <- function(df_tp) {
    df_tp %>%
        group_by(subject, feature) %>%
        summarise(value = value[which.max(!is.na(value))], .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = value) %>%
        arrange(subject) %>%
        as.data.frame() -> wide
    rownames(wide) <- wide$subject
    wide[, setdiff(colnames(wide), "subject"), drop = FALSE]
}

dat_long <- prep_long(combined_exeromedat_long)

# Enforce allowed times per fluid explicitly (belt & suspenders)
allowed_times <- list(
    plasma = c("-1","0","0.5","1","3","24"),
    saliva = c("-1","0","0.5","1","3","24"),
    urine  = c("-1","0","24")
)
dat_long <- dat_long %>%
    filter(mapply(function(f,t) t %in% allowed_times[[f]], fluid, time))

# Timepoints present in data (and in desired order)
all_tp <- intersect(names(time_color_palette), unique(dat_long$time))
timepoints <- factor(all_tp, levels = names(time_color_palette)) |> as.character()
message("Timepoints found: ", paste(timepoints, collapse = ", "))

# Output dirs and clean stale outputs
unlink(list.files(FILES_DIR,
                  pattern="^(filtered_correlations_|directionally_unique_correlations_|cor_full_matrix_|pval_full_matrix_|nobs_matrix_|widemat_timepoint_).*(rds|csv|png)$",
                  full.names=TRUE))
unlink(list.files(FIG_DIR, pattern="^(network_plot_|fluid_pair_counts|hist_cor_).*(pdf|png)$", full.names=TRUE))

# We'll collect raw pairwise correlations across timepoints for global histograms
hist_source_all <- list()

# Global histogram fallback thresholds — define ONCE up here
PAIR_N_MIN_HIST <- if (exists("PAIR_N_MIN_HIST")) PAIR_N_MIN_HIST else 3

# ------------------------ Pairwise Spearman with guards ---------------------
for (tp in timepoints) {
    cat("Processing timepoint:", tp, "\n")

    # Filter to the fluids allowed at this tp
    allowed_fluids <- allowed_fluids_by_tp[[tp]]
    tp_df <- dat_long %>% filter(time == tp, fluid %in% allowed_fluids)
    if (nrow(tp_df) == 0) next

    mat <- tp_to_matrix(tp_df)

    # Coerce to numeric and ensure clean colnames
    mat[] <- lapply(mat, function(v) suppressWarnings(as.numeric(v)))
    mat <- as.data.frame(mat, check.names = TRUE)

    # Assert no list columns and all numeric
    if (any(vapply(mat, is.list, logical(1)))) {
        stop("List-columns detected after pivot; expected plain numeric columns.")
    }
    if (!all(vapply(mat, is.numeric, logical(1)))) {
        bad <- names(mat)[!vapply(mat, is.numeric, logical(1))]
        stop("Non-numeric columns remain after coercion: ", paste(bad, collapse = ", "))
    }

    # Drop columns with < N_MIN finite values OR no variance
    valid_cols <- sapply(mat, function(v) {
        v <- v[is.finite(v)]
        length(v) >= N_MIN && length(unique(v)) >= 2
    })
    if (!all(valid_cols)) mat <- mat[, valid_cols, drop = FALSE]
    if (ncol(mat) < 2) {
        warning("Not enough valid features at timepoint ", tp, " for correlations.")
        next
    }

    # Build containers
    n <- ncol(mat)
    cor_matrix <- matrix(NA_real_, n, n, dimnames = list(colnames(mat), colnames(mat)))
    p_matrix   <- matrix(NA_real_, n, n, dimnames = list(colnames(mat), colnames(mat)))
    n_matrix   <- matrix(NA_integer_, n, n, dimnames = list(colnames(mat), colnames(mat)))

    # ------- PARALLEL pairwise correlations (upper triangle) -------
    pairs_ij <- utils::combn(n, 2L)
    if (PARALLEL_COR) {
        res <- foreach(k = seq_len(ncol(pairs_ij)),
                       .combine = rbind,
                       .inorder = FALSE) %dopar% {
                           i <- pairs_ij[1, k]; j <- pairs_ij[2, k]
                           x <- mat[[i]]; y <- mat[[j]]
                           ok   <- is.finite(x) & is.finite(y)
                           n_ok <- sum(ok)
                           if (n_ok < N_MIN ||
                               length(unique(x[ok])) < 2 ||
                               length(unique(y[ok])) < 2) {
                               data.frame(i=i, j=j, r=NA_real_, p=NA_real_, n=n_ok)
                           } else {
                               ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
                               data.frame(i=i, j=j, r=unname(ct$estimate), p=ct$p.value, n=n_ok)
                           }
                       }
    } else {
        res <- do.call(rbind, lapply(seq_len(ncol(pairs_ij)), function(k) {
            i <- pairs_ij[1, k]; j <- pairs_ij[2, k]
            x <- mat[[i]]; y <- mat[[j]]
            ok   <- is.finite(x) & is.finite(y)
            n_ok <- sum(ok)
            if (n_ok < N_MIN ||
                length(unique(x[ok])) < 2 ||
                length(unique(y[ok])) < 2) {
                data.frame(i=i, j=j, r=NA_real_, p=NA_real_, n=n_ok)
            } else {
                ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
                data.frame(i=i, j=j, r=unname(ct$estimate), p=ct$p.value, n=n_ok)
            }
        }))
    }
    if (!is.null(res) && nrow(res)) {
        idx <- cbind(res$i, res$j)
        cor_matrix[idx] <- res$r; cor_matrix[cbind(res$j, res$i)] <- res$r
        p_matrix[idx]   <- res$p; p_matrix[cbind(res$j, res$i)]   <- res$p
        n_matrix[idx]   <- res$n; n_matrix[cbind(res$j, res$i)]   <- res$n
    }
    diag(cor_matrix) <- 1
    diag(p_matrix)   <- 1
    diag(n_matrix)   <- sapply(mat, function(v) sum(is.finite(v)))

    # Adjust p-values (upper triangle)
    upper_idx <- which(upper.tri(p_matrix))
    adj_pvals <- p.adjust(p_matrix[upper_idx], method = PADJ_METHOD)
    adj_p_matrix <- matrix(NA_real_, n, n, dimnames = list(colnames(mat), colnames(mat)))
    adj_p_matrix[upper_idx] <- adj_pvals
    adj_p_matrix[lower.tri(adj_p_matrix)] <- t(adj_p_matrix)[lower.tri(adj_p_matrix)]
    diag(adj_p_matrix) <- 1

    # ---- QA: verify correlation estimates on a random subset (optional) ----
    if (CHECK_COR) {
        idxu <- which(upper.tri(cor_matrix) & is.finite(cor_matrix))
        if (length(idxu) > 0) {
            sample_idx <- sample(idxu, size = min(CHECK_COR_N, length(idxu)))
            diffs <- numeric(length(sample_idx))
            for (k in seq_along(sample_idx)) {
                ij <- arrayInd(sample_idx[k], .dim = c(n, n))
                i <- ij[1]; j <- ij[2]
                x <- mat[[i]]; y <- mat[[j]]
                ok <- is.finite(x) & is.finite(y)
                if (sum(ok) >= N_MIN) {
                    r1 <- suppressWarnings(cor(x[ok], y[ok], method="spearman"))
                    r2 <- cor_matrix[i, j]
                    diffs[k] <- abs(r1 - r2)
                }
            }
            max_diff <- max(diffs, na.rm = TRUE)
            if (is.finite(max_diff) && max_diff > 1e-6) {
                warning(sprintf("QA: Max |Δρ| between cor() and cor.test() at tp %s = %.3g", tp, max_diff))
            } else {
                message(sprintf("QA: Spearman estimates consistent at tp %s (max |Δρ| ≤ 1e-6).", tp))
            }
        }
    }

    # Save core matrices + wide mat
    saveRDS(cor_matrix,   file = p_files(sprintf("cor_full_matrix_timepoint_%s.rds", tp)))
    saveRDS(adj_p_matrix, file = p_files(sprintf("pval_full_matrix_timepoint_%s_adjusted.rds", tp)))
    saveRDS(n_matrix,     file = p_files(sprintf("nobs_matrix_timepoint_%s.rds", tp)))
    saveRDS(mat,          file = p_files(sprintf("widemat_timepoint_%s.rds", tp)))

    # ----- Build long pairs_df and per-timepoint histograms (padj<0.05 in red) -----
    ut <- which(upper.tri(cor_matrix), arr.ind = TRUE)
    pairs_df <- data.frame(
        Molecule1 = rownames(cor_matrix)[ut[,1]],
        Molecule2 = colnames(cor_matrix)[ut[,2]],
        R         = cor_matrix[ut],
        P_adj     = adj_p_matrix[ut],
        N_pair    = n_matrix[ut],
        stringsAsFactors = FALSE
    ) %>%
        dplyr::filter(is.finite(R), N_pair >= N_MIN) %>%
        dplyr::mutate(
            Fluid1 = get_fluid_from_feature(Molecule1),
            Fluid2 = get_fluid_from_feature(Molecule2)
        ) %>%
        dplyr::filter(Fluid1 %in% allowed_fluids, Fluid2 %in% allowed_fluids) %>%
        dplyr::mutate(
            FluidPair = paste(pmin(Fluid1, Fluid2), pmax(Fluid1, Fluid2), sep = FLUID_SEP),
            timepoint = tp
        ) %>%
        dplyr::filter(!CROSS_FLUID_ONLY | Fluid1 != Fluid2)

    if (nrow(pairs_df)) {
        hist_source_all[[tp]] <- pairs_df

        sig_counts <- pairs_df %>%
            dplyr::group_by(FluidPair) %>%
            dplyr::summarise(n_sig = sum(P_adj < 0.05, na.rm = TRUE), .groups = "drop")

        p_hist <- ggplot(pairs_df, aes(x = R)) +
            geom_histogram(binwidth = 0.05, boundary = 0, closed = "right") +
            geom_histogram(
                data = subset(pairs_df, P_adj < 0.05),
                binwidth = 0.05, boundary = 0, closed = "right",
                fill = "red", alpha = 0.45
            ) +
            geom_text(
                data = sig_counts,
                aes(x = Inf, y = Inf, label = paste0("n(p_adj<0.05) = ", n_sig)),
                inherit.aes = FALSE, color = "red", size = 3.2,
                hjust = 1.1, vjust = 1.2
            ) +
            geom_vline(xintercept = c(-R_MIN, R_MIN), linetype = 2) +
            facet_wrap(~ FluidPair, scales = "free_y") +
            labs(title = paste("Distribution of Spearman rho by Fluid Pair -- tp", tp),
                 x = "Spearman rho", y = "Count") +
            theme_minimal()

        ggsave(filename = p_fig(sprintf("hist_cor_by_fluidpair_tp_%s.pdf", tp)),
               plot = p_hist, width = 8, height = 5,
               device = if (USE_CAIRO && "cairo_pdf" %in% ls("package:grDevices"))
                   grDevices::cairo_pdf else "pdf")
    }

    # ----- Global histograms across all timepoints (if any) -----
    # Upper triangle indices
    ut <- which(upper.tri(cor_matrix), arr.ind = TRUE)

    pairs_df <- data.frame(
        Molecule1 = rownames(cor_matrix)[ut[,1]],
        Molecule2 = colnames(cor_matrix)[ut[,2]],
        R         = cor_matrix[ut],
        P_adj     = adj_p_matrix[ut],
        N_pair    = n_matrix[ut],
        stringsAsFactors = FALSE
    ) %>%
        dplyr::filter(is.finite(R), N_pair >= PAIR_N_MIN_HIST) %>%   # <-- looser pairwise N
        dplyr::mutate(
            Fluid1 = get_fluid_from_feature(Molecule1),
            Fluid2 = get_fluid_from_feature(Molecule2)
        ) %>%
        dplyr::filter(Fluid1 %in% allowed_fluids, Fluid2 %in% allowed_fluids) %>%
        dplyr::mutate(
            FluidPair = paste(pmin(Fluid1, Fluid2), pmax(Fluid1, Fluid2), sep = FLUID_SEP),
            timepoint = tp
        ) %>%
        dplyr::filter(!CROSS_FLUID_ONLY | Fluid1 != Fluid2)

    if (nrow(pairs_df)) {
        # Ensure ASCII-only facet labels (avoid en-dashes/emdashes)
        pairs_df$FluidPair <- gsub("[\u2012-\u2015\u2212]", FLUID_SEP, pairs_df$FluidPair)

        # Keep for global histograms
        if (!exists("hist_source_all")) hist_source_all <- list()
        hist_source_all[[tp]] <- pairs_df

        # Count padj<0.05 per facet
        sig_counts <- pairs_df %>%
            dplyr::group_by(FluidPair) %>%
            dplyr::summarise(n_sig = sum(P_adj < 0.05, na.rm = TRUE), .groups = "drop")

        # Plot: base histogram + red overlay for padj<0.05 + red count label
        p_hist <- ggplot(pairs_df, aes(x = R)) +
            geom_histogram(binwidth = 0.05, boundary = 0, closed = "right") +
            geom_histogram(
                data = subset(pairs_df, P_adj < 0.05),
                binwidth = 0.05, boundary = 0, closed = "right",
                fill = "red", alpha = 0.45
            ) +
            geom_text(
                data = sig_counts,
                aes(x = Inf, y = Inf, label = paste0("n(p_adj<0.05) = ", n_sig)),
                inherit.aes = FALSE, color = "red", size = 3.2,
                hjust = 1.1, vjust = 1.2
            ) +
            geom_vline(xintercept = c(-R_MIN, R_MIN), linetype = 2) +
            facet_wrap(~ FluidPair, scales = "free_y") +
            labs(title = paste("Distribution of Spearman rho by Fluid Pair -- tp", tp),
                 x = "Spearman rho", y = "Count") +
            theme_minimal()

        # Save PDF
        out_path <- p_fig(sprintf("hist_cor_by_fluidpair_tp_%s.pdf", tp))
        ggsave(filename = out_path, plot = p_hist, width = 8, height = 5,
               device = if (USE_CAIRO && "cairo_pdf" %in% ls("package:grDevices"))
                   grDevices::cairo_pdf else "pdf")
        message("Wrote: ", normalizePath(out_path, winslash = "/"))

        message(sprintf("tp %s: histogram pairs kept = %d (>= %d overlap)",
                        tp, nrow(pairs_df), PAIR_N_MIN_HIST))
    } else {
        message(sprintf("tp %s: no pairs >= %d overlap; skipping histogram.",
                        tp, PAIR_N_MIN_HIST))
    }

}  # <-- CLOSES the outer for (tp in timepoints) { ... } loop

# ==================== EDGE LISTS + BOOTSTRAP STABILITY ======================
for (tp in timepoints) {
    cat("Annotating edges for:", tp, "\n")

    cor_matrix   <- readRDS(p_files(sprintf("cor_full_matrix_timepoint_%s.rds", tp)))
    adj_p_matrix <- readRDS(p_files(sprintf("pval_full_matrix_timepoint_%s_adjusted.rds", tp)))
    n_matrix     <- readRDS(p_files(sprintf("nobs_matrix_timepoint_%s.rds", tp)))
    widemat      <- readRDS(p_files(sprintf("widemat_timepoint_%s.rds", tp)))

    if (is.null(dim(cor_matrix)) || is.null(dim(adj_p_matrix))) next

    cor_df <- as.data.frame(as.table(cor_matrix))
    names(cor_df) <- c("Molecule1", "Molecule2", "Correlation")

    # Drop diagonal + duplicate undirected pairs
    cor_df <- cor_df %>% filter(Molecule1 != Molecule2)
    cor_df <- cor_df[!duplicated(apply(cor_df[,1:2], 1, function(x) paste(sort(x), collapse = ""))), ]

    # Bring adjusted p-value and pairwise N
    cor_df$adj_p_value <- adj_p_matrix[cbind(cor_df$Molecule1, cor_df$Molecule2)]
    cor_df$N_pair      <- n_matrix[cbind(cor_df$Molecule1, cor_df$Molecule2)]

    # Thresholds for edges
    filtered_cor <- cor_df %>%
        filter(
            is.finite(Correlation),
            is.finite(adj_p_value),
            N_pair >= N_MIN,
            adj_p_value < ALPHA_ADJ,
            abs(Correlation) >= R_MIN
        ) %>%
        arrange(desc(abs(Correlation))) %>%
        mutate(
            Fluid1 = get_fluid_from_feature(Molecule1),
            Fluid2 = get_fluid_from_feature(Molecule2)
        ) %>%
        filter(!CROSS_FLUID_ONLY | Fluid1 != Fluid2)

    # ------- PARALLEL bootstrap stability; keep only stable edges -------
    if (nrow(filtered_cor)) {
        if (PARALLEL_BOOT) {
            filtered_cor$BootStability <- furrr::future_pmap_dbl(
                list(filtered_cor$Molecule1, filtered_cor$Molecule2),
                function(m1, m2) {
                    x <- widemat[[m1]]; y <- widemat[[m2]]
                    boot_stability_pair(x, y, B = B_BOOT, r_min = R_MIN, n_min = N_MIN)
                },
                .options = furrr::furrr_options(seed = TRUE)
            )
        } else {
            filtered_cor$BootStability <- purrr::pmap_dbl(
                list(filtered_cor$Molecule1, filtered_cor$Molecule2),
                function(m1, m2) {
                    x <- widemat[[m1]]; y <- widemat[[m2]]
                    boot_stability_pair(x, y, B = B_BOOT, r_min = R_MIN, n_min = N_MIN)
                }
            )
        }
        filtered_cor <- filtered_cor %>% filter(BootStability >= STABILITY_MIN)
    }

    # Derive fluids from names and enforce allowed for this tp (defensive)
    if (nrow(filtered_cor)) {
        allowed <- allowed_fluids_by_tp[[tp]]
        bad <- which(!(filtered_cor$Fluid1 %in% allowed & filtered_cor$Fluid2 %in% allowed))
        if (length(bad)) {
            message("Dropping ", length(bad), " edges with disallowed fluids at tp ", tp)
            filtered_cor <- filtered_cor[-bad, , drop = FALSE]
        }
        if (tp %in% c("0.5","1","3") &&
            any(filtered_cor$Fluid1 == "urine" | filtered_cor$Fluid2 == "urine")) {
            offenders <- head(filtered_cor[filtered_cor$Fluid1 == "urine" | filtered_cor$Fluid2 == "urine",
                                           c("Molecule1","Molecule2","Correlation")], 5)
            stop("Sanity check failed: urine edges present at tp ", tp, ". Examples:\n",
                 paste(capture.output(print(offenders)), collapse = "\n"))
        }
    }

    write.csv(filtered_cor,
              file = p_files(sprintf("filtered_correlations_%s.csv", tp)),
              row.names = FALSE)
}

# =================== DIRECTIONALLY UNIQUE FILTER ============================
cor_files <- list.files(FILES_DIR, pattern = "^filtered_correlations_.*\\.csv$", full.names = TRUE)
cor_list <- lapply(cor_files, function(p) read.csv(p, check.names = FALSE))
names(cor_list) <- gsub("filtered_correlations_|\\.csv", "", basename(cor_files))

make_pair_key <- function(m1, m2) apply(cbind(m1, m2), 1, function(x) paste(sort(x), collapse = "_"))

for (tp in names(cor_list)) {
    df <- cor_list[[tp]]
    if (!nrow(df)) next
    df$Fluid1 <- get_fluid_from_feature(df$Molecule1)
    df$Fluid2 <- get_fluid_from_feature(df$Molecule2)
    allowed <- allowed_fluids_by_tp[[tp]]
    df <- df[df$Fluid1 %in% allowed & df$Fluid2 %in% allowed, , drop = FALSE]
    if (tp %in% c("0.5","1","3") && nrow(df) &&
        any(df$Fluid1 == "urine" | df$Fluid2 == "urine")) {
        stop("Sanity check failed in directionally-unique stage at tp ", tp, ".")
    }
    df$pair_key <- make_pair_key(df$Molecule1, df$Molecule2)
    cor_list[[tp]] <- df
}

all_pairs_long <- do.call(rbind, lapply(names(cor_list), function(tp) {
    df <- cor_list[[tp]]
    if (is.null(df) || !nrow(df)) return(NULL)
    data.frame(
        timepoint     = tp,
        pair_key      = df$pair_key,
        Correlation   = df$Correlation,
        Molecule1     = df$Molecule1,
        Molecule2     = df$Molecule2,
        adj_p_value   = df$adj_p_value,
        BootStability = if ("BootStability" %in% names(df)) df$BootStability else NA_real_,
        Fluid1        = df$Fluid1,
        Fluid2        = df$Fluid2,
        stringsAsFactors = FALSE
    )
}))

if (!is.null(all_pairs_long) && nrow(all_pairs_long) > 0) {
    all_pairs_long$sign <- sign(all_pairs_long$Correlation)

    pair_signs  <- aggregate(sign ~ pair_key, data = all_pairs_long, function(x) length(unique(x)))
    pair_counts <- aggregate(timepoint ~ pair_key, data = all_pairs_long, function(x) length(unique(x)))

    pair_summary <- merge(pair_signs, pair_counts, by = "pair_key")
    names(pair_summary) <- c("pair_key", "n_signs", "n_timepoints")
    pair_summary$keep <- (pair_summary$n_timepoints == 1) | (pair_summary$n_signs > 1)

    valid_pairs <- pair_summary$pair_key[pair_summary$keep]
    final_filtered_df <- all_pairs_long[all_pairs_long$pair_key %in% valid_pairs, ]

    split_filtered <- split(final_filtered_df, final_filtered_df$timepoint)
    for (tp in names(split_filtered)) {
        df <- split_filtered[[tp]]
        if (tp %in% c("0.5","1","3") && any(df$Fluid1 == "urine" | df$Fluid2 == "urine")) {
            stop("Sanity check failed just before writing directionally_unique_correlations_", tp, ".csv")
        }
        write.csv(df,
                  file = p_files(sprintf("directionally_unique_correlations_%s.csv", tp)),
                  row.names = FALSE)
    }
}

# ============================== NETWORKS ====================================
all_network_metrics <- data.frame()
for (tp in timepoints) {
    cat("Network for:", tp, "\n")
    f <- p_files(sprintf("directionally_unique_correlations_%s.csv", tp))
    if (!file.exists(f)) next
    df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    if (!nrow(df)) next

    # Recompute fluids and enforce allowed (defensive)
    df$Fluid1 <- get_fluid_from_feature(df$Molecule1)
    df$Fluid2 <- get_fluid_from_feature(df$Molecule2)
    allowed <- allowed_fluids_by_tp[[tp]]
    df <- df[df$Fluid1 %in% allowed & df$Fluid2 %in% allowed, , drop = FALSE]
    if (tp %in% c("0.5","1","3") && any(df$Fluid1 == "urine" | df$Fluid2 == "urine")) {
        stop("Sanity check failed in NETWORKS stage at tp ", tp, ".")
    }

    needed <- c("Molecule1","Molecule2","Correlation")
    if (!all(needed %in% names(df))) {
        stop("Missing columns in ", f, ": ", paste(setdiff(needed, names(df)), collapse = ", "))
    }

    edge_df <- df %>%
        dplyr::select(Molecule1, Molecule2,
                      Correlation,
                      dplyr::any_of(c("BootStability","adj_p_value","Fluid1","Fluid2")))

    g <- graph_from_data_frame(edge_df, directed = FALSE)
    g <- igraph::simplify(
        g,
        remove.multiple = TRUE,
        remove.loops    = TRUE,
        edge.attr.comb  = list(
            Correlation   = "mean",
            BootStability = "mean",
            adj_p_value   = "min",
            weight        = "mean",
            stability     = "mean",
            .default      = "first"
        )
    )

    E(g)$weight    <- E(g)$Correlation
    E(g)$stability <- if ("BootStability" %in% names(edge_df)) E(g)$BootStability else 1
    E(g)$color     <- ifelse(E(g)$weight > 0, "#ADD8E6", "#ef8a62")

    base_w <- 0.5 + 5.5 * rescale01(abs(E(g)$weight))   # ~0.5–6
    E(g)$width <- base_w * (0.5 + 0.5 * E(g)$stability)

    all_nodes <- unique(c(edge_df$Molecule1, edge_df$Molecule2))
    node_fluids <- setNames(c(df$Fluid1, df$Fluid2),
                            c(df$Molecule1, df$Molecule2))[all_nodes]
    V(g)$color <- fluid_color_palette[node_fluids[V(g)$name]]

    deg <- degree(g)
    hub_cutoff <- stats::quantile(deg, 0.80, na.rm = TRUE)
    is_hub <- deg > hub_cutoff

    V(g)$shape <- ifelse(is_hub[V(g)$name], "square", "circle")
    V(g)$size  <- 3 + 6 * rescale01(deg)        # ~3–9
    V(g)$size[is_hub[V(g)$name]] <- V(g)$size[is_hub[V(g)$name]] + 1.5
    V(g)$label <- NA
    V(g)$label.cex <- 0.5

    set.seed(123)
    l <- layout_with_fr(g, weights = abs(E(g)$weight))

    pdf(p_fig(sprintf("network_plot_%s.pdf", tp)), width = 12/2.54, height = 9/2.54)
    plot(g,
         layout = l,
         main = paste("Cross-Fluid Protein Network at", tp),
         vertex.label.color = "black",
         vertex.frame.color = NA,
         edge.curved = FALSE,
         margin = c(0,0,0,0))
    legend("bottomright", legend = names(fluid_color_palette), col = fluid_color_palette, pch = 16,
           title = "Fluid", pt.cex = 0.7, bty = "n", cex = 0.7)
    legend("right", legend = c("Hub", "Other"), pch = c(15,16), col = "black",
           title = "Node", pt.cex = 0.7, bty = "n", cex = 0.7)
    legend("bottomleft", legend = c("Positive", "Negative"), lty = 1,
           col = c("#ADD8E6","#ef8a62"), lwd = 2, bty = "n", cex = 0.7)
    mtext(paste("Nodes:", vcount(g), "| Edges:", ecount(g),
                "| Hubs:", sum(is_hub)), side = 1, line = 0, cex = 0.7)
    dev.off()

    all_network_metrics <- rbind(all_network_metrics, data.frame(
        Timepoint = tp,
        Nodes = vcount(g),
        Edges = ecount(g),
        Hubs = sum(is_hub),
        Density = edge_density(g),
        AvgDegree = mean(degree(g))
    ))
}
write.csv(all_network_metrics, p_files("network_metrics.csv"), row.names = FALSE)

# ======================== FLUID-PAIR COUNTS PLOT ============================
pair_counts <- data.frame()
for (tp in timepoints) {
    f <- p_files(sprintf("directionally_unique_correlations_%s.csv", tp))
    if (!file.exists(f)) next
    df <- read.csv(f, check.names = FALSE)
    if (!nrow(df)) next

    df$Fluid1 <- get_fluid_from_feature(df$Molecule1)
    df$Fluid2 <- get_fluid_from_feature(df$Molecule2)

    allowed <- allowed_fluids_by_tp[[tp]]
    bad <- which(!(df$Fluid1 %in% allowed & df$Fluid2 %in% allowed))
    if (length(bad)) {
        message("PairCounts: dropping ", length(bad), " disallowed rows at tp ", tp)
        df <- df[-bad, , drop = FALSE]
    }
    if (tp %in% c("0.5","1","3") && any(df$Fluid1 == "urine" | df$Fluid2 == "urine")) {
        stop("Sanity check failed in pair-counts at tp ", tp, ".")
    }

    if (!nrow(df)) next
    df$FluidPair <- apply(df[, c("Fluid1","Fluid2")], 1, function(x) paste(sort(x), collapse = FLUID_SEP))
    counts <- as.data.frame(table(df$FluidPair))
    colnames(counts) <- c("FluidPair", "Count")
    counts$Timepoint <- tp
    pair_counts <- rbind(pair_counts, counts)
}
if (nrow(pair_counts) > 0) {
    total_counts <- aggregate(Count ~ FluidPair + Timepoint, data = pair_counts, sum)
    pdf(p_fig("fluid_pair_counts.pdf"), width = 6, height = 4)
    print(
        ggplot(total_counts, aes(x = reorder(FluidPair, -Count), y = Count)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            facet_wrap(~ Timepoint, scales = "free_y") +
            labs(title = "Counts of Fluid Pairs by Timepoint",
                 x = "Fluid Pair", y = "Count") +
            theme_minimal() +
            coord_flip()
    )
    dev.off()
}


# ============================ PART 2 — build networks from the saved matrices ==
# Reads the saved correlation / p-value / n-obs matrices, filters cross-fluid
# edges (padj < ALPHA_ADJ, N_pair >= N_MIN; no |r| threshold), detects Louvain
# modules, runs per-module GO:BP ORA (a CSV is written only for significant
# modules), and plots modules of size >= MIN_MODULE_SIZE_FOR_PLOT.

suppressPackageStartupMessages({   # clusterProfiler / org.Hs.eg.db for the ORA (others loaded at top)
    library(clusterProfiler); library(org.Hs.eg.db)
})

set.seed(1234)  # re-seed so Part 2 (Louvain, ORA) is reproducible independently of Part 1

# --------------------------- Config / paths ---------------------------
if (!exists("PROJECT_ROOT")) PROJECT_ROOT <- normalizePath(".", mustWork = FALSE)
OUT_DIR   <- if (exists("OUT_DIR")) OUT_DIR else PROJECT_ROOT
FILES_DIR <- if (exists("FILES_DIR")) FILES_DIR else file.path(OUT_DIR, "files")
FIG_DIR   <- if (exists("FIG_DIR"))   FIG_DIR   else file.path(OUT_DIR, "figures")
if (!dir.exists(FILES_DIR)) dir.create(FILES_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(FIG_DIR))   dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)
p_files <- function(...) file.path(FILES_DIR, ...)
p_fig   <- function(...) file.path(FIG_DIR, ...)

ALPHA_ADJ <- if (exists("ALPHA_ADJ")) ALPHA_ADJ else 0.05
N_MIN     <- if (exists("N_MIN"))     N_MIN     else 16

fluid_color_palette <- if (exists("fluid_color_palette")) fluid_color_palette else c(
    "plasma"="#4B0082","saliva"="#FFA07A","urine"="#FFD700"
)
allowed_fluids_by_tp <- if (exists("allowed_fluids_by_tp")) allowed_fluids_by_tp else list(
    "-1"=c("plasma","saliva","urine"), "0"=c("plasma","saliva","urine"),
    "0.5"=c("plasma","saliva"), "1"=c("plasma","saliva"), "3"=c("plasma","saliva"),
    "24"=c("plasma","saliva","urine")
)

# timepoints from saved matrices (unless provided)
if (!exists("timepoints")) {
    tp_files <- list.files(FILES_DIR, pattern="^cor_full_matrix_timepoint_.*\\.rds$", full.names=FALSE)
    timepoints <- sort(unique(sub("^cor_full_matrix_timepoint_(.*)\\.rds$", "\\1", tp_files)))
}

# ------------------------------ Helpers --------------------------------
get_fluid_from_feature <- function(s) { f <- sub(".*_", "", s); tolower(sub("\\..*$","", trimws(f))) }
strip_fluid_suffix <- function(x) sub("_(plasma|saliva|urine)$", "", x, ignore.case = TRUE)
node_is_prot <- function(n) grepl("_(plasma|saliva|urine)$", n, ignore.case = TRUE)
rescale01 <- function(x){ xr <- range(x, na.rm=TRUE); if(!all(is.finite(xr))||diff(xr)==0) return(rep(0.5,length(x))); (x-xr[1])/diff(xr) }

# Build cross-fluid, padj-only edge table from saved matrices
build_cf_edges_for_tp <- function(tp, alpha=ALPHA_ADJ, n_min=N_MIN){
    cor_matrix   <- readRDS(p_files(sprintf("cor_full_matrix_timepoint_%s.rds", tp)))
    adj_p_matrix <- readRDS(p_files(sprintf("pval_full_matrix_timepoint_%s_adjusted.rds", tp)))
    n_matrix     <- readRDS(p_files(sprintf("nobs_matrix_timepoint_%s.rds", tp)))
    ut <- which(upper.tri(cor_matrix), arr.ind = TRUE)
    df <- data.frame(
        Molecule1   = rownames(cor_matrix)[ut[,1]],
        Molecule2   = colnames(cor_matrix)[ut[,2]],
        Correlation = as.numeric(cor_matrix[ut]),
        adj_p_value = as.numeric(adj_p_matrix[ut]),
        N_pair      = as.integer(n_matrix[ut]),
        stringsAsFactors = FALSE
    )
    df$Fluid1 <- get_fluid_from_feature(df$Molecule1)
    df$Fluid2 <- get_fluid_from_feature(df$Molecule2)
    allowed <- allowed_fluids_by_tp[[tp]]
    df <- df[df$Fluid1 %in% allowed & df$Fluid2 %in% allowed, , drop=FALSE]
    df <- df[
        df$Fluid1 != df$Fluid2 &
            is.finite(df$adj_p_value) & df$adj_p_value < alpha &
            is.finite(df$N_pair) & df$N_pair >= n_min &
            is.finite(df$Correlation),
        , drop=FALSE
    ]
    df
}

# One canonical membership CSV per timepoint
cf_membership_path <- function(tp) p_files(sprintf("module_membership_crossfluid_padj_%s.csv", tp))
get_or_make_cf_membership <- function(g, tp){
    f <- cf_membership_path(tp)
    if (file.exists(f)){
        mem_df <- tryCatch(read.csv(f, stringsAsFactors=FALSE), error=function(e) NULL)
        if (!is.null(mem_df) && all(c("Molecule","Module") %in% names(mem_df))) {
            mm <- setNames(as.integer(mem_df$Module), mem_df$Molecule)
            memb <- mm[V(g)$name]
            if (!any(is.na(memb))) return(memb)
        }
    }
    set.seed(123)
    lou  <- cluster_louvain(g, weights = abs(E(g)$weight))
    memb <- as.integer(membership(lou)); names(memb) <- V(g)$name
    write.csv(data.frame(Molecule = names(memb), Module = memb, Timepoint = tp),
              f, row.names = FALSE)
    memb
}

# ORA universe (simple version, as you had it)
get_assays <- function(prot_label_df){
    ac <- names(prot_label_df)[match(tolower("Assay"), tolower(names(prot_label_df)))]
    if (is.na(ac)) stop("prot.label is missing an 'Assay' column.")
    as.character(trimws(prot_label_df[[ac]]))
}
build_universe_for_tp <- function(tp){
    fluids <- if (tp %in% c("-1","0","24")) c("plasma","saliva","urine") else c("plasma","saliva")
    assays <- character(0)
    if ("plasma" %in% fluids && exists("olink.prot.label"))   assays <- c(assays, get_assays(olink.prot.label)) # these can be any list or gene names from any file
    if ("saliva" %in% fluids && exists("prot.label.saliva"))  assays <- c(assays, get_assays(prot.label.saliva))
    if ("urine"  %in% fluids && exists("prot.label.urine"))   assays <- c(assays, get_assays(prot.label.urine))
    unique(assays[nzchar(assays)])
}


# ------------------------------------------------------------------
# DISTINCT MODULE PALETTE - complex in case there are many modules
# ------------------------------------------------------------------
module_palette <- function(k, seed = 123,
                           # HCL bands to sample from (tweak if needed):
                           L_vals = c(70, 80),     # lightness (higher = lighter)
                           C_vals = c(55, 70),     # chroma (higher = more vivid)
                           hue_step = 4            # candidate hue granularity (degrees)
) {
    if (k <= 0) return(character(0))
    set.seed(seed)

    # 1) Build candidate pool in HCL (hue across the wheel, 2×L, 2×C bands)
    hues <- seq(0, 356, by = hue_step)  # 0..356
    cand <- expand.grid(h = hues, L = L_vals, C = C_vals, KEEP.OUT.ATTRS = FALSE)
    # Remove invalid (out-of-gamut) colors
    cand$hex <- grDevices::hcl(h = cand$h, c = cand$C, l = cand$L)
    is_ok <- !is.na(cand$hex)
    cand  <- cand[is_ok, , drop = FALSE]

    # If still too few candidates, widen bands
    if (nrow(cand) < k) {
        L_vals2 <- sort(unique(c(L_vals, mean(range(L_vals)) + c(-10, 10))))
        C_vals2 <- sort(unique(c(C_vals, mean(range(C_vals)) + c(-10, 10))))
        hues2   <- seq(0, 355, by = max(2, hue_step %/% 2))
        cand2 <- expand.grid(h = hues2, L = L_vals2, C = C_vals2, KEEP.OUT.ATTRS = FALSE)
        cand2$hex <- grDevices::hcl(h = cand2$h, c = cand2$C, l = cand2$L)
        cand <- unique(rbind(cand, cand2[!is.na(cand2$hex), , drop = FALSE]))
    }

    # 2) Greedy farthest-point sampling in HCL (Euclidean in HCL works well)
    #    Start from a random candidate with highest chroma (vivid seed)
    idx_pool <- seq_len(nrow(cand))
    pool_vivid <- idx_pool[order(-cand$C, -cand$L)][seq_len(min(20, nrow(cand)))]
    sel_idx <- sample(pool_vivid, 1)
    selected <- cand[sel_idx, c("h","C","L"), drop = FALSE]

    # helper: distance from each candidate to nearest selected color
    dist_to_sel <- function(M, S) {
        # M: matrix of candidates (h,C,L); S: matrix of selected (h,C,L)
        # circular hue distance: wrap via min(|dh|, 360-|dh|)
        dh <- abs(outer(M[,1], S[,1], FUN = "-"))
        dh <- pmin(dh, 360 - dh)
        dC <- abs(outer(M[,2], S[,2], FUN = "-"))
        dL <- abs(outer(M[,3], S[,3], FUN = "-"))
        # weight hue a bit higher; tune weights if you like
        w_h <- 1.2; w_c <- 1.0; w_l <- 0.8
        D <- sqrt((w_h*dh)^2 + (w_c*dC)^2 + (w_l*dL)^2)
        apply(D, 1, min)
    }

    remaining <- setdiff(idx_pool, sel_idx)
    while (length(sel_idx) < min(k, nrow(cand))) {
        M <- as.matrix(cand[remaining, c("h","C","L")])
        S <- as.matrix(selected)
        dmin <- dist_to_sel(M, S)
        pick <- remaining[ which.max(dmin) ]
        sel_idx <- c(sel_idx, pick)
        selected <- rbind(selected, cand[pick, c("h","C","L")])
        remaining <- setdiff(remaining, pick)
    }

    cand$hex[sel_idx][seq_len(k)]
}

# ----------------- change module size (number of proteins) ----------------

EDGE_POS_COL <- "#0072B2"; EDGE_NEG_COL <- "#D55E00"
MIN_MODULE_SIZE_FOR_ORA  <- 10L
MIN_MODULE_SIZE_FOR_PLOT <- 10L

# ===================== Loop A: ORA (only save if significant) =====================
FORCE_ORA <- TRUE  # set FALSE to keep existing significant files

for (tp in timepoints) {
    cat("ORA pass for tp:", tp, "\n")
    df <- tryCatch(build_cf_edges_for_tp(tp), error=function(e){message("  Missing matrices for tp ", tp); NULL})
    if (is.null(df) || !nrow(df)) { message("  No significant cross-fluid pairs; skipping."); next }

    edge_df <- df %>%
        dplyr::select(Molecule1, Molecule2, Correlation, Fluid1, Fluid2) %>%
        dplyr::filter(Molecule1 != Molecule2)
    edge_df <- edge_df[!duplicated(t(apply(edge_df[,1:2], 1, sort))), , drop=FALSE]
    if (!nrow(edge_df)) { message("  No edges to draw; skipping."); next }

    g <- graph_from_data_frame(edge_df, directed = FALSE)
    E(g)$weight <- as.numeric(E(g)$Correlation); E(g)$weight[!is.finite(E(g)$weight)] <- 0

    memb <- get_or_make_cf_membership(g, tp)
    V(g)$Module <- as.integer(memb)
    tab <- table(V(g)$Module)
    keep_mods <- as.integer(names(tab)[tab >= MIN_MODULE_SIZE_FOR_ORA])
    if (!length(keep_mods)) { message("  No modules ≥ ", MIN_MODULE_SIZE_FOR_ORA, " for ORA."); next }

    universe <- build_universe_for_tp(tp)

    for (mid in keep_mods) {
        out_sig <- p_files(sprintf("ORA_sig_crossfluid_padj_%s_Module_%s.csv", tp, mid))
        if (!FORCE_ORA && file.exists(out_sig)) next

        mem_nodes  <- V(g)$name[V(g)$Module == mid]
        prot_nodes <- mem_nodes[node_is_prot(mem_nodes)]
        prot_syms  <- unique(strip_fluid_suffix(prot_nodes))
        genes      <- prot_syms[prot_syms %in% universe]

        if (length(genes) < 3) { if (file.exists(out_sig)) file.remove(out_sig); next }

        ora <- tryCatch({
            clusterProfiler::enrichGO(
                gene          = genes,
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2
            )
        }, error=function(e) NULL)

        if (is.null(ora)) { if (file.exists(out_sig)) file.remove(out_sig); next }

        df_ora <- as.data.frame(ora)
        if (!nrow(df_ora)) { if (file.exists(out_sig)) file.remove(out_sig); next }

        # keep only significant rows (by ALPHA_ADJ) to decide saving
        lc <- tolower(names(df_ora))
        pcol <- match(TRUE, lc %in% c("p.adjust","padj","p_adj","qvalue","q.value"))
        if (is.na(pcol)) { if (file.exists(out_sig)) file.remove(out_sig); next }
        pv <- suppressWarnings(as.numeric(df_ora[[pcol]]))
        sig_rows <- is.finite(pv) & pv < ALPHA_ADJ
        if (any(sig_rows)) {
            write.csv(df_ora[sig_rows, , drop=FALSE], out_sig, row.names=FALSE)  # <<< ONLY save when significant
        } else if (file.exists(out_sig)) {
            file.remove(out_sig)  # ensure no stale file
        }
    }
}

# ===================== Loop B: Plot (uses ORA_sig files for coloring) =====================
all_network_metrics_cf <- data.frame()

for (tp in timepoints) {
    cat("Plotting tp:", tp, "\n")
    df <- tryCatch(build_cf_edges_for_tp(tp), error=function(e){message("  Missing matrices for tp ", tp); NULL})
    if (is.null(df) || !nrow(df)) { message("  No significant cross-fluid pairs; skipping."); next }

    edge_df <- df %>%
        dplyr::select(Molecule1, Molecule2, Correlation, Fluid1, Fluid2) %>%
        dplyr::filter(Molecule1 != Molecule2)
    edge_df <- edge_df[!duplicated(t(apply(edge_df[,1:2], 1, sort))), , drop=FALSE]
    if (!nrow(edge_df)) { message("  No edges to draw; skipping."); next }

    g <- graph_from_data_frame(edge_df, directed = FALSE)
    E(g)$weight <- as.numeric(E(g)$Correlation); E(g)$weight[!is.finite(E(g)$weight)] <- 0
    E(g)$color  <- ifelse(E(g)$weight > 0, EDGE_POS_COL, EDGE_NEG_COL)
    E(g)$width  <- 0.5 + 5.5 * rescale01(abs(E(g)$weight))

    # node attrs
    nodes <- unique(c(edge_df$Molecule1, edge_df$Molecule2))
    node_fluids <- setNames(c(df$Fluid1, df$Fluid2), c(df$Molecule1, df$Molecule2))[nodes]
    V(g)$fluid <- node_fluids[V(g)$name]
    V(g)$color <- fluid_color_palette[V(g)$fluid]
    V(g)$frame.color <- "#333333"; V(g)$frame.width <- 0.5

    # membership
    memb <- get_or_make_cf_membership(g, tp)
    V(g)$Module <- as.integer(memb)

    # keep only modules with >= MIN_MODULE_SIZE_FOR_PLOT
    tab <- table(V(g)$Module)
    plot_mods <- as.integer(names(tab)[tab >= MIN_MODULE_SIZE_FOR_PLOT])
    if (!length(plot_mods)) { message("  No modules ≥ ", MIN_MODULE_SIZE_FOR_PLOT, "; skipping plot."); next }
    g_plot <- induced_subgraph(g, vids = which(V(g)$Module %in% plot_mods))

    # hubs/labels on g_plot
    deg_plot  <- degree(g_plot)
    strg_plot <- strength(g_plot, vids = V(g_plot), weights = abs(E(g_plot)$weight))
    hub_cut   <- stats::quantile(deg_plot, 0.80, na.rm=TRUE)
    is_hub    <- deg_plot > hub_cut
    V(g_plot)$shape <- ifelse(is_hub[V(g_plot)$name], "square", "circle")
    V(g_plot)$size  <- 3 + 6 * rescale01(deg_plot)
    V(g_plot)$size[is_hub[V(g_plot)$name]] <- V(g_plot)$size[is_hub[V(g_plot)$name]] + 1.5
    V(g_plot)$label <- NA; V(g_plot)$label.cex <- 0.65
    pick_labels <- function(gp, k=5L){
        labs <- character(); d <- degree(gp); s <- strength(gp, weights=abs(E(gp)$weight))
        for (fl in unique(V(gp)$fluid)) {
            vv <- V(gp)$name[V(gp)$fluid==fl]; if (!length(vv)) next
            ord <- order(-d[vv], -s[vv], vv); labs <- c(labs, head(vv[ord], k))
        }
        unique(labs)
    }
    to_lab <- pick_labels(g_plot, 5L)
    V(g_plot)$label[match(to_lab, V(g_plot)$name)] <- strip_fluid_suffix(to_lab)

    # module groups and ORA significance mask from presence of ORA_sig file
    kept_mods <- sort(unique(V(g_plot)$Module))
    groups <- lapply(kept_mods, function(m) which(V(g_plot)$Module == m))
    names(groups) <- paste0("M", kept_mods)

    ora_sig_flags <- logical(length(kept_mods)); names(ora_sig_flags) <- kept_mods
    for (i in seq_along(kept_mods)) {
        f <- p_files(sprintf("ORA_sig_crossfluid_padj_%s_Module_%s.csv", tp, kept_mods[i]))
        ora_sig_flags[i] <- file.exists(f)
    }

    # colors: ORA+ use distinct Set3; ORA- light grey
    k_sig <- sum(ora_sig_flags)
    pal_sig <- module_palette(k_sig)
    mark_cols   <- character(length(kept_mods))
    mark_border <- character(length(kept_mods))
    sig_idx <- which(ora_sig_flags)
    if (k_sig > 0) {
        mark_cols[sig_idx]   <- grDevices::adjustcolor(pal_sig, alpha.f = 0.40)
        mark_border[sig_idx] <- pal_sig
    }
    non_idx <- which(!ora_sig_flags)
    if (length(non_idx)) {
        mark_cols[non_idx]   <- grDevices::adjustcolor("#B0B0B0", alpha.f = 0.18)
        mark_border[non_idx] <- NA
    }

    # layout/plot
    set.seed(123)
    lay <- layout_with_fr(g_plot, weights = abs(E(g_plot)$weight))

    pdf(p_fig(sprintf("network_plot_crossfluid_padj_ORAaware_%s.pdf", tp)), width = 12/2.54, height = 9/2.54)
    par(mar = c(2, 0.5, 2, 0.5))
    plot(g_plot,
         layout = lay,
         vertex.color       = V(g_plot)$color,
         vertex.frame.color = V(g_plot)$frame.color,
         vertex.shape       = V(g_plot)$shape,
         vertex.size        = V(g_plot)$size,
         vertex.label       = V(g_plot)$label,
         vertex.label.color = "black",
         vertex.label.cex   = V(g_plot)$label.cex,
         edge.width         = E(g_plot)$width,
         edge.color         = ifelse(E(g_plot)$weight>0, EDGE_POS_COL, EDGE_NEG_COL),
         edge.curved        = FALSE,
         mark.groups        = groups,
         mark.col           = mark_cols,
         mark.border        = mark_border,
         mark.expand        = 4,
         margin             = c(0,0,0,0))
    title(main = paste0("Cross-Fluid Network at ", tp), cex.main = 0.8, line = 0.2)

    # Legend only for ORA+ modules
    if (k_sig > 0) {
        legend("topright",
               legend = paste0("M", kept_mods[ora_sig_flags]),
               pch = 15, pt.cex = 1, bty = "n", cex = 0.7,
               col = mark_cols[ora_sig_flags], title = "Sig. ORA modules")
    }
    legend("bottomleft", legend=c("Positive","Negative"), lty=1,
           col=c(EDGE_POS_COL, EDGE_NEG_COL), lwd=2, bty="n", cex=0.7, title="Edge sign")
    legend("bottomright", legend=names(fluid_color_palette), col=fluid_color_palette,
           pch=16, title="Body Fluid", pt.cex=0.7, bty="n", cex=0.7)

    mtext(paste("Nodes:", vcount(g_plot), "| Edges:", ecount(g_plot), "| Hubs:", sum(is_hub)),
          side = 1, line = 0, cex = 0.7)
    dev.off()

    # exports/metrics
    igraph::write_graph(g_plot, p_files(sprintf("network_crossfluid_padj_ORAaware_%s.graphml", tp)), format = "graphml")
    write.csv(data.frame(
        Molecule = V(g_plot)$name,
        LabelStripped = strip_fluid_suffix(V(g_plot)$name),
        Module  = V(g_plot)$Module,
        Fluid   = V(g_plot)$fluid,
        Degree  = as.numeric(deg_plot[V(g_plot)]),
        Strength= as.numeric(strg_plot[V(g_plot)]),
        IsHub   = is_hub[V(g_plot)]
    ), p_files(sprintf("network_crossfluid_padj_ORAaware_node_attributes_%s.csv", tp)), row.names = FALSE)

    all_network_metrics_cf <- rbind(all_network_metrics_cf, data.frame(
        Timepoint = tp,
        Nodes = vcount(g_plot),
        Edges = ecount(g_plot),
        Hubs = sum(is_hub),
        Density = edge_density(g_plot),
        AvgDegree = mean(deg_plot)
    ))
}
write.csv(all_network_metrics_cf, p_files("network_metrics_crossfluid_padj_ORAaware.csv"), row.names = FALSE)

