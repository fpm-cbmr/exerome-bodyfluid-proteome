# =============================================================================
# Extended Data Fig. 1 (related to main Fig. 4) — Cross-fluid Spearman correlation networks at Pre and 0.5 h, 3 h, and 24 h post-exercise
#
# Input:  data/network/network_crossfluid_padj_ORAaware_<tp>.graphml
#           <- 14_crossfluid_correlation_networks.R
# Output: figures/extended_data_1/*.pdf
# =============================================================================
suppressMessages({ library(here); library(igraph); library(writexl) })
source(here("R/figure_defaults.R"))
source(here("R/figures/helpers_crossfluid_network.R"))
NET_DIR <- here("data", "network")
FIG <- here("figures", "extended_data_1"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(here("source_data"), showWarnings = FALSE)

panels <- list(a = "-1", b = "0.5", c = "3", d = "24")
tp_title <- c(`-1` = "Pre", `0.5` = "0.5 h", `3` = "3 h", `24` = "24 h")
SHEETS <- list()
for (pn in names(panels)) {
  tp <- panels[[pn]]
  g <- cf_plot_network(tp, file.path(FIG, sprintf("%s_network_%s.pdf", pn, tp)), NET_DIR,
                       title = paste0("Cross-Fluid Network at ", tp_title[tp]))
  SHEETS[[sprintf("%s_network_%s_inset", pn, tp)]] <-
    cf_top_hub_inset(g, file.path(FIG, sprintf("%s_network_%s_inset.pdf", pn, tp)))
  SHEETS[[sprintf("%s_network_%s_edges", pn, tp)]] <- igraph::as_data_frame(g, "edges")
}
# source data now built by R/source_data/_source_data.R
cat("Extended Data 1 panels written:", paste(names(SHEETS), collapse = ", "), "\n")
