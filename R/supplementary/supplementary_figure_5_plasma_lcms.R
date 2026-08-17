# =============================================================================
# Supplementary Figure 5 — Plasma proteome response to acute exercise measured with LC-MS/MS proteomics (related to Fig. 3).
#
#   a-d  temporally resolved LC-MS/MS plasma clusters; lowest-BH-adj-P proteins
#        coloured, all others grey  [function = cluster_profile_lcms(), baseline-normalized data]
#   e    Olink vs LC-MS correlation for the 170 shared proteins vs mean LC-MS
#        abundance, exercise-regulated proteins highlighted red(+)/blue(-)
#        (produced by supplementary_figure_5e_olink_lcms_exercise_highlight.R -> S5e_*)
#
# Input:  data/res.plasma.linear.rda, data/exerome.dat.plasma.rda (baseline-
#         normalized LC-MS), data/prot.label.plasma.rda
# Output: figures/supplementary/S5_*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2); library(ggpubr); library(ggsci)
})
source(here("R", "package_loading.R")); source(here("R", "functions_loading.R"))
OUT <- here("figures", "supplementary"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

load(here("data", "res.plasma.linear.rda"))
load(here("data", "exerome.dat.plasma.rda")); exerome.dat <- exerome.dat.plasma   # baseline-normalized LC-MS
load(here("data", "prot.label.plasma.rda"));  prot.label <- prot.label.plasma

cl.label.exerome <- tibble::tibble(
  cluster = 1:4,
  label = c("Increase, over correction at 24h", "Decrease to 24h",
            "Elevated in recovery", "Stochastic patterns"),
  cluster.paper = paste0("C", 1:4))

source(here("R", "supplementary", "helpers_supp_plots.R"))   # write_supp_source_data()
SHEETS <- list()
# ---- a-d: per-cluster temporal profiles ------------
for (k in 1:4) {
  # LC-MS z-scores span roughly -3..2
  pk <- cluster_profile_lcms(res.plasma.linear, exerome.dat, prot.label, cl.label.exerome, cluster_number = k, ylim = c(-3, 2))
  ggsave(file.path(OUT, sprintf("S5%s_plasma_lcms_cluster_%d.pdf", letters[k], k)), pk,
         width = 5.5, height = 5, units = "cm", dpi = 600)
  SHEETS[[sprintf("%s_lcms_cluster%d", letters[k], k)]] <- pk$data
}
write_supp_source_data(SHEETS, "5")   # merges with panel e (supplementary_figure_5e_olink_lcms_exercise_highlight.R)
cat("Supplementary Figure 5 a-d (LC-MS clusters) written; panel e from supplementary_figure_5e_olink_lcms_exercise_highlight.R\n")
