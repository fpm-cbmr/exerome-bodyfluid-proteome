# =============================================================================
# Plasma Olink — linear mixed model + temporal clusters + variance decomposition
#                (Discovery cohort).
#
# Baseline (t = -1) z-normalisation of the Olink NPX, then a linear mixed model
# per assay (random subject intercept; age + VO2max covariates). Assays with a
# time effect at fdr.aov < 0.05 are grouped into temporal clusters with navmix
# (K = 5, seed 42).
# Also computes, per protein, the between-subject baseline variance (t = -1) and
# the within-subject residual variance (Fig 3b).
#
# Input:  data/olink_npx.data.rda, data/olink.prot.label.rda, data-raw/clinical_data.xlsx
# Output: data/res.olink.linear.rda,
#         data/olink.exerome.dat.rda, data/plasma_variance_decomposition.rda
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr)
  library(data.table); library(doParallel); library(parallel); library(lme4); library(navmix)
})
source(here("R/package_loading.R")); source(here("R/functions_loading.R"))

load(here("data", "olink_npx.data.rda"))
npx.data <- olink_npx.data %>% dplyr::rename(subject = label)      # cols: OID*, time, label
load(here("data", "olink.prot.label.rda")); prot.label <- olink.prot.label
clinical_data <- readxl::read_xlsx(here("data-raw", "clinical_data.xlsx"))
clinical_data$subject <- gsub("^EX", "", clinical_data$subject)

# baseline z-normalisation
exerome.norm <- npx.data %>% dplyr::filter(time == -1) %>%
  tidyr::pivot_longer(cols = starts_with("OID"), names_to = "OlinkID", values_to = "value") %>%
  dplyr::group_by(OlinkID) %>%
  dplyr::summarize(mean.value = mean(value, na.rm = TRUE), sd.value = sd(value, na.rm = TRUE), .groups = "drop")
exerome.dat <- npx.data
for (j in seq_len(nrow(exerome.norm)))
  exerome.dat <- exerome.dat %>%
    dplyr::mutate(!!exerome.norm$OlinkID[j] := (!!sym(exerome.norm$OlinkID[j]) - exerome.norm$mean.value[j]) / exerome.norm$sd.value[j])
exerome.dat <- exerome.dat %>% left_join(clinical_data, by = "subject") %>% dplyr::mutate(t.factor = factor(time))
olink.exerome.dat <- exerome.dat %>% dplyr::mutate(time_label = as.numeric(as.character(time))) %>% dplyr::select(-time)
save(olink.exerome.dat, file = here("data", "olink.exerome.dat.rda"))

# linear mixed model
registerDoParallel(cores = max(1, detectCores() - 1))
res.olink.linear <- mixed_anova_parallel(exerome.dat, "t.factor", unique(prot.label$OlinkID),
                                         "+ (1|subject)", covariates = c("age", "vo2_max"))
stopImplicitCluster()
res.olink.linear <- res.olink.linear %>%
  left_join(prot.label, by = c("outcome" = "OlinkID")) %>%
  mutate(fdr.aov = p.adjust(pval.aov.t.factor, method = "BH"))
res.olink.linear <- as.data.table(res.olink.linear) %>% dplyr::rename(OlinkID = outcome) %>%
  dplyr::distinct(OlinkID, .keep_all = TRUE)

# navmix temporal clustering (relabelled to the published clusters below)
cl.data <- res.olink.linear %>% dplyr::filter(fdr.aov < .05) %>%
  dplyr::select(OlinkID, matches("beta|se|pval")) %>% as.data.frame()
rownames(cl.data) <- cl.data$OlinkID
tps <- c(0, 0.5, 1, 3, 24)
for (j in tps) cl.data[, paste0("zscore.", j)] <- cl.data[, paste0("beta.exposure", j)] / cl.data[, paste0("se.exposure", j)]
set.seed(42)
fit <- navmix(as.matrix(dplyr::select(cl.data, starts_with("zscore."))), K = 5, plot = FALSE, plot_radial = FALSE)
cl.data$cluster <- fit$fit$z
# navmix labels are arbitrary, so give each recomputed cluster the number of the
# published cluster (Supplementary Table 2, sheet D) that its proteins most overlap.
pub_cl <- as.data.frame(suppressWarnings(readxl::read_xlsx(
  here("data-raw", "supplementary_table_2.xlsx"), sheet = "D")))
pub_cl <- setNames(suppressWarnings(as.integer(pub_cl$cluster)), pub_cl$OlinkID)
ov <- table(cl.data$cluster, pub_cl[rownames(cl.data)])
relabel <- setNames(as.integer(colnames(ov)[max.col(ov, ties.method = "first")]), rownames(ov))
tmp.cl <- data.frame(OlinkID = rownames(cl.data), cluster = unname(relabel[as.character(cl.data$cluster)]))
res.olink.linear <- res.olink.linear %>% left_join(tmp.cl, by = "OlinkID") %>% arrange(fdr.aov)
save(res.olink.linear, file = here("data", "res.olink.linear.rda"))
sig <- res.olink.linear %>% filter(fdr.aov < 0.05)
cat(sprintf("Plasma: %d assays, %d exercise-regulated (%d unique genes). Cluster sizes:\n",
            nrow(res.olink.linear), nrow(sig), dplyr::n_distinct(sig$Assay)))
print(table(res.olink.linear$cluster))

# ---- variance decomposition (Fig 3b): baseline vs residual variance ----------
protein_columns <- intersect(prot.label$OlinkID, grep("^OID", names(npx.data), value = TRUE))
npx_long <- npx.data %>% pivot_longer(cols = all_of(protein_columns), names_to = "Assay", values_to = "value")
npx_filtered <- npx_long %>% group_by(Assay, subject) %>% filter(n() >= 3) %>% ungroup()
calc_baseline <- function(a, data) {
  d <- subset(data, Assay == a & time == -1)
  data.frame(Assay = a, Baseline_Variance = if (nrow(d) < 2) NA else var(d$value, na.rm = TRUE))
}
calc_residual <- function(a, data) {
  d <- subset(data, Assay == a); if (nrow(d) < 2) return(data.frame(Assay = a, Residual_Variance = NA))
  m <- tryCatch(lme4::lmer(value ~ time + (1 | subject), data = d, REML = TRUE), error = function(e) NULL)
  data.frame(Assay = a, Residual_Variance = if (is.null(m)) NA else as.data.frame(lme4::VarCorr(m))[2, "vcov"])
}
cl <- makeCluster(max(1, detectCores() - 1))
clusterExport(cl, c("npx_filtered", "calc_baseline", "calc_residual"), envir = environment())
clusterEvalQ(cl, library(lme4))
bv <- do.call(rbind, parLapply(cl, unique(npx_filtered$Assay), calc_baseline, data = npx_filtered))
rv <- do.call(rbind, parLapply(cl, unique(npx_filtered$Assay), calc_residual, data = npx_filtered))
stopCluster(cl)
plasma_variance_decomposition <- merge(bv, rv, by = "Assay") %>% dplyr::rename(OlinkID = Assay) %>%
  left_join(res.olink.linear, by = "OlinkID")
save(plasma_variance_decomposition, file = here("data", "plasma_variance_decomposition.rda"))
cat("Plasma variance decomposition saved (Fig 3b).\n")
