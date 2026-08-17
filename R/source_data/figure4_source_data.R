# =============================================================================
# Source Data — Figure 4 (cross-fluid networks). Built from committed data/ analysis
# objects: the cross-fluid Spearman networks (data/network/*.graphml, analysis 14)
# and the cross-fluid container (data/se_bodyfluid.rda). One sheet per panel;
# network panels export the node + edge tables straight from the graphml.
#
# Output: source_data/SourceData_Figure4.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(stringr)
  library(igraph); library(SummarizedExperiment); library(writexl) })
source(here("R/figures/helpers_crossfluid_network.R"))   # cf_graphml_path
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()
NET_DIR <- here("data", "network")
tps <- c("-1", "0", "0.5", "1", "3", "24"); tp_labels <- c("Pre", "0h", "0.5h", "1h", "3h", "24h")

# ---- b: total vs cross-fluid edge counts per timepoint -----------------------
metrics <- read.csv(file.path(NET_DIR, "network_metrics.csv"))
cf_edges <- sapply(tps, function(tp) { f <- cf_graphml_path(tp, NET_DIR)
  if (file.exists(f)) igraph::ecount(igraph::read_graph(f, "graphml")) else NA })
SHEETS[["b_edge_counts"]] <- data.frame(
  Timepoint = tp_labels,
  TotalEdges = metrics$Edges[match(as.numeric(tps), metrics$Timepoint)],
  CrossFluidEdges = as.integer(cf_edges))

# ---- c, f: cross-fluid protein correlations from se_bodyfluid ----------------
load(here("data", "se_bodyfluid.rda"))
se_data <- as.data.frame(SummarizedExperiment::assay(se_bodyfluid))
filtered_data <- se_data[, grepl("Pre|T0|T24", colnames(se_data))]
core_names <- gsub("_(urine|saliva|plasma)", "", rownames(filtered_data))
matching_cores <- names(table(core_names)[table(core_names) > 1])
filtered_data <- filtered_data[rownames(filtered_data)[core_names %in% matching_cores], ]

mansc1 <- filtered_data[grep("MANSC1", rownames(filtered_data)), ]
SHEETS[["c_MANSC1"]] <- data.frame(
  Saliva = as.numeric(mansc1[grepl("_saliva", rownames(mansc1)), ]),
  Plasma = as.numeric(mansc1[grepl("_plasma", rownames(mansc1)), ]),
  Time_Point = sub(".*_", "", colnames(mansc1))) %>% dplyr::filter(!is.na(Saliva) & !is.na(Plasma))

protein_names <- rownames(filtered_data)
res_rows <- list()
for (core in matching_cores) {
  pr <- filtered_data[grep(paste0("^", core, "_"), protein_names), , drop = FALSE]
  if (nrow(pr) < 2) next
  for (tm in c("Pre", "T0", "T24")) {
    tcols <- grepl(tm, colnames(filtered_data)); if (!sum(tcols)) next
    td <- pr[, tcols, drop = FALSE]
    for (i in 1:(nrow(td) - 1)) for (j in (i + 1):nrow(td)) {
      x <- as.numeric(td[i, ]); y <- as.numeric(td[j, ])
      if (sum(!is.na(x) & !is.na(y)) < 2) next
      tt <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
      src <- paste0(sub(".*_", "", rownames(pr)[i]), "_", sub(".*_", "", rownames(pr)[j]))
      res_rows[[length(res_rows) + 1]] <- data.frame(Protein = core, Correlation = unname(tt$estimate),
                                                     Source = src, Time_Point = tm)
    }
  }
}
SHEETS[["f_rho_pre_vs_0h"]] <- dplyr::bind_rows(res_rows) %>%
  dplyr::mutate(Time_Point = dplyr::recode(Time_Point, T0 = "0h", T24 = "24h")) %>%
  dplyr::filter(!Source %in% c("urine_urine", "saliva_saliva")) %>%
  dplyr::mutate(Source = stringr::str_replace(Source, "olink", "plasma")) %>%
  dplyr::group_by(Protein, Source, Time_Point) %>%
  dplyr::summarise(Correlation = mean(Correlation, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(Time_Point %in% c("Pre", "0h")) %>%
  tidyr::pivot_wider(id_cols = c(Protein, Source), names_from = Time_Point,
                     values_from = Correlation, names_prefix = "Time_") %>%
  tidyr::drop_na(Time_Pre, Time_0h) %>% as.data.frame()

# ---- d, e: cross-fluid network node + edge tables (0 h, 1 h) ------------------
for (tp in c("0", "1")) {
  lab <- if (tp == "0") "d_network_0h" else "e_network_1h"
  g <- igraph::read_graph(cf_graphml_path(tp, NET_DIR), "graphml")
  SHEETS[[paste0(lab, "_nodes")]] <- igraph::as_data_frame(g, "vertices")
  SHEETS[[paste0(lab, "_edges")]] <- igraph::as_data_frame(g, "edges")
}

writexl::write_xlsx(SHEETS, here("source_data", "SourceData_Figure4.xlsx"))
cat("Source Data Figure 4 written:", paste(names(SHEETS), collapse = ", "), "\n")
