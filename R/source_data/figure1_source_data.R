# =============================================================================
# Source Data — Figure 1.
#
# Builds SourceData_Figure1.xlsx directly from the committed analysis objects in
# data/ (produced by R/analysis/), independent of the figure-rendering script.
# Each sheet is the table underlying one panel.
#
# Input:  data/res.msd.linear.rda, data/msd.exerome.dat.rda,
#         data/metabolic.cytokine.prot.label.rda, data/res.saliva.linear.rda,
#         data/res.urine.linear.rda, data/res.olink.linear.rda
# Output: source_data/SourceData_Figure1.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)

load(here("data", "res.msd.linear.rda"))
load(here("data", "msd.exerome.dat.rda"))
load(here("data", "res.saliva.linear.rda"))
load(here("data", "res.urine.linear.rda"))
load(here("data", "res.olink.linear.rda"))

# ---- b, c: per-timepoint mean z-score for the metabolic (1) / cytokine (2) sets
exerome.dat <- msd.exerome.dat %>% dplyr::rename(time_label = time)
marker_profiles <- function(cluster_num) {
  mts <- res.msd.linear %>% dplyr::filter(cluster == cluster_num) %>% dplyr::pull(Assay)
  dplyr::bind_rows(lapply(mts, function(m) {
    exerome.dat %>% dplyr::group_by(time_label) %>%
      dplyr::summarise(mean_zscore = mean(.data[[m]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(Assay = m)
  })) %>% dplyr::mutate(cluster = cluster_num) %>%
    dplyr::select(Assay, cluster, time_label, mean_zscore)
}

# ---- d: significant-protein counts per fluid / timepoint / direction ---------
count_sig <- function(df, fluid_name) {
  df <- as.data.frame(df)
  pv <- grep("^pval\\.exposure", names(df), value = TRUE)
  for (p in pv) df[[paste0(p, "_fdr")]] <- p.adjust(df[[p]], "fdr")
  fdr <- df %>% dplyr::select(Assay, dplyr::matches("pval\\.exposure.*_fdr")) %>%
    tidyr::pivot_longer(-Assay, names_to = "time_point", names_pattern = "pval\\.exposure(.*)_fdr", values_to = "pval_fdr")
  bet <- df %>% dplyr::select(Assay, dplyr::matches("^beta\\.exposure")) %>%
    tidyr::pivot_longer(-Assay, names_to = "time_point", names_pattern = "beta\\.exposure(.*)", values_to = "beta")
  fdr %>% dplyr::inner_join(bet, by = c("Assay", "time_point"), relationship = "many-to-many") %>%
    dplyr::filter(pval_fdr < 0.05) %>%
    dplyr::mutate(direction = ifelse(beta > 0, "increasing", "decreasing"), fluid = fluid_name) %>%
    dplyr::group_by(fluid, time_point, direction) %>% dplyr::summarise(count = dplyr::n(), .groups = "drop")
}
panel_d <- dplyr::bind_rows(count_sig(res.saliva.linear, "Saliva"),
                            count_sig(res.urine.linear,  "Urine"),
                            count_sig(res.olink.linear,  "Plasma")) %>%
  dplyr::mutate(time_point = factor(time_point, levels = c("0", "0.5", "1", "3", "24")))

# ---- e: proteins shared across >=2 fluids, by directional concordance --------
fluid_direction <- function(res, fluid_name) {
  res <- as.data.frame(res)
  beta_cols <- grep("^beta\\.exposure", colnames(res), value = TRUE)
  sig <- res[res$fdr.aov < 0.05 & !is.na(res$fdr.aov), ]
  peak <- apply(as.matrix(sig[, beta_cols, drop = FALSE]), 1, function(b) {
    b <- b[!is.na(b)]; if (length(b) == 0) NA_real_ else b[which.max(abs(b))] })
  tibble::tibble(Assay = sig$Assay, fluid = fluid_name, direction = ifelse(peak > 0, "up", "down")) %>%
    dplyr::filter(!is.na(direction)) %>% dplyr::distinct(Assay, .keep_all = TRUE)
}
panel_e <- dplyr::bind_rows(fluid_direction(res.olink.linear, "Plasma"),
                            fluid_direction(res.saliva.linear, "Saliva"),
                            fluid_direction(res.urine.linear, "Urine")) %>%
  dplyr::group_by(Assay) %>%
  dplyr::summarise(fluids = paste(sort(fluid), collapse = "+"), n_fluids = dplyr::n(),
                   n_up = sum(direction == "up"), n_down = sum(direction == "down"), .groups = "drop") %>%
  dplyr::filter(n_fluids >= 2) %>%
  dplyr::mutate(Direction = dplyr::case_when(n_down == 0 ~ "Concordant up",
                                             n_up == 0   ~ "Concordant down",
                                             TRUE        ~ "Discordant"))

lineage <- data.frame(
  panel = c("b", "c", "d", "e"),
  description = c("Plasma metabolic markers over time", "MSD cytokine responses over time",
                 "No. significantly regulated proteins per fluid/timepoint",
                 "Proteins shared across >=2 fluids by directional concordance (UpSet)"),
  analysis_script = c("01_metabolic_cytokine_lmm.R", "01_metabolic_cytokine_lmm.R",
                      "02/03/03b/04 *_lmm_clusters.R", "02/03/04 *_lmm_clusters.R"),
  source_object = c("res.msd.linear + msd.exerome.dat", "res.msd.linear + msd.exerome.dat",
                    "res.{saliva,urine,olink}.linear", "res.{olink,saliva,urine}.linear"))

writexl::write_xlsx(list(lineage = lineage,
                         b_metabolic_markers = marker_profiles(1),
                         c_cytokines         = marker_profiles(2),
                         d_protein_counts    = as.data.frame(panel_d),
                         e_upset_shared      = as.data.frame(panel_e)),
                    here("source_data", "SourceData_Figure1.xlsx"))
cat("Source Data Figure 1 written to source_data/SourceData_Figure1.xlsx\n")
