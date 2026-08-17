# =============================================================================
# Source Data — Extended Data Fig. 1 (cross-fluid networks at Pre/0.5/3/24 h).
# Node + edge tables exported straight from the committed network graphml objects
# (data/network/, analysis 14). One node+edge sheet pair per panel/timepoint.
#
# Output: source_data/SourceData_ExtendedData1.xlsx
# =============================================================================
suppressMessages({ library(here); library(igraph); library(writexl) })
source(here("R/figures/helpers_crossfluid_network.R"))   # cf_graphml_path
dir.create(here("source_data"), showWarnings = FALSE)
NET_DIR <- here("data", "network")
panels <- list(a = "-1", b = "0.5", c = "3", d = "24")
SHEETS <- list()
for (pn in names(panels)) {
  tp <- panels[[pn]]
  g <- igraph::read_graph(cf_graphml_path(tp, NET_DIR), "graphml")
  SHEETS[[sprintf("%s_network_%s_nodes", pn, tp)]] <- igraph::as_data_frame(g, "vertices")
  SHEETS[[sprintf("%s_network_%s_edges", pn, tp)]] <- igraph::as_data_frame(g, "edges")
}
writexl::write_xlsx(SHEETS, here("source_data", "SourceData_ExtendedData1.xlsx"))
cat("Source Data Extended Data 1 written:", paste(names(SHEETS), collapse = ", "), "\n")
