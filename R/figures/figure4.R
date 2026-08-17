# =============================================================================
# Figure 4 — Acute exercise coordinates protein networks across body fluids.
#
# Panels:
#   a    study schematic (BioRender; no source data)
#   b    total vs cross-fluid correlation-edge counts over time
#   c    MANSC1 cross-fluid example (saliva vs plasma z-score, per time point)
#   d,e  cross-fluid correlation networks at 0 h and 1 h post-exercise (+ top-hub insets)
#   f    per-protein cross-fluid Spearman correlations, pre- vs 0 h (3 fluid pairs)
#
# Inputs (produced by the analysis pipeline):
#   data/network/network_metrics.csv                            <- 14_crossfluid_correlation_networks.R   (b)
#   data/network/network_crossfluid_padj_ORAaware_<tp>.graphml  <- 14_crossfluid_correlation_networks.R   (b,d,e)
#   data/se_bodyfluid.rda                                       <- 10_bodyfluid_data_container.R          (c,f)
#
# Output: figures/figure_4/*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggpubr); library(ggrepel); library(patchwork)
  library(writexl); library(igraph); library(SummarizedExperiment)
})
source(here("R/package_loading.R")); source(here("R/functions_loading.R")); source(here("R/figure_defaults.R"))
source(here("R/figures/helpers_crossfluid_network.R"))
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
count <- dplyr::count; rename <- dplyr::rename; summarise <- dplyr::summarise
group_by <- dplyr::group_by; arrange <- dplyr::arrange; slice <- dplyr::slice
dir.create(here("figures", "figure_4"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

NET_DIR      <- here("data", "network")
tps          <- c("-1", "0", "0.5", "1", "3", "24")
tp_labels    <- c("Pre", "0h", "0.5h", "1h", "3h", "24h")
fluid_color_palette <- c(plasma = "#4B0082", saliva = "#FFA07A", urine = "#FFD700")
EDGE_POS_COL <- "#0072B2"; EDGE_NEG_COL <- "#D55E00"

graphml_path <- function(tp) cf_graphml_path(tp, NET_DIR)

# ---- panel b: total vs cross-fluid edges per timepoint -----------------------
metrics <- read.csv(file.path(NET_DIR, "network_metrics.csv"))
cf_edges <- sapply(tps, function(tp) {
  f <- graphml_path(tp); if (file.exists(f)) ecount(read_graph(f, "graphml")) else NA
})
edge_df <- data.frame(Timepoint = factor(tp_labels, levels = tp_labels),
                      TotalEdges = metrics$Edges[match(as.numeric(tps), metrics$Timepoint)],
                      CrossFluidEdges = as.integer(cf_edges))
sc <- max(edge_df$TotalEdges, na.rm = TRUE) / max(edge_df$CrossFluidEdges, na.rm = TRUE)
pb <- ggplot(edge_df, aes(x = Timepoint, group = 1)) +
  geom_line(aes(y = TotalEdges)) + geom_point(aes(y = TotalEdges), shape = 16, size = 1.4) +
  geom_line(aes(y = CrossFluidEdges * sc), linetype = "dashed") +
  geom_point(aes(y = CrossFluidEdges * sc), shape = 21, fill = "white", size = 1.4) +
  scale_y_continuous(name = "Total edges",
                     sec.axis = sec_axis(~ . / sc, name = "Cross-fluid edges")) +
  labs(x = "Timepoint", title = "protein-protein correlations") +
  theme_pubr(base_size = 6) + theme(text = element_text(size = 6)) + theme_strokes
ggsave(here("figures", "figure_4", "b_edge_counts.pdf"), pb, width = 6, height = 5, units = "cm", dpi = 600)
SHEETS[["b_edge_counts"]] <- edge_df

# =============================================================================
# panels c & f — cross-fluid protein correlations from the WGCNA SE
# =============================================================================
load(here("data", "se_bodyfluid.rda"))
se_data       <- as.data.frame(assay(se_bodyfluid))
filtered_data <- se_data[, grepl("Pre|T0|T24", colnames(se_data))]
protein_names <- rownames(filtered_data)
core_names    <- gsub("_(urine|saliva|plasma)", "", protein_names)
matching_cores <- names(table(core_names)[table(core_names) > 1])
filtered_data <- filtered_data[protein_names[core_names %in% matching_cores], ]
protein_names <- rownames(filtered_data)

# ---- panel c: MANSC1 saliva vs plasma z-score, per time point ----------------
mansc1 <- filtered_data[grep("MANSC1", rownames(filtered_data)), ]
c_df <- data.frame(
  Saliva = as.numeric(mansc1[grepl("_saliva", rownames(mansc1)), ]),
  Plasma = as.numeric(mansc1[grepl("_plasma", rownames(mansc1)), ]),
  Time_Point = sub(".*_", "", colnames(mansc1))
) %>% filter(!is.na(Saliva) & !is.na(Plasma))
c_time_pal <- c(Pre = "#E69F00", `T0` = "#56B4E9", `T24` = "#CC79A7")
c_stats <- c_df %>% group_by(Time_Point) %>%
  summarise(Correlation = cor(Saliva, Plasma, method = "spearman"), .groups = "drop")
pc <- ggplot(c_df, aes(Saliva, Plasma)) +
  geom_point(aes(color = Time_Point), size = 1, alpha = 0.7) +
  geom_smooth(aes(group = Time_Point, color = Time_Point, fill = Time_Point),
              method = "lm", se = TRUE, size = 0.5, alpha = 0.3) +
  scale_color_manual(values = c_time_pal) + scale_fill_manual(values = c_time_pal) +
  geom_text(data = c_stats, inherit.aes = FALSE, size = 2, hjust = 0.5,
            aes(x = mean(c_df$Saliva, na.rm = TRUE) + c(-1, 0, 1)[match(Time_Point, c("Pre","T0","T24"))],
                y = max(c_df$Plasma, na.rm = TRUE) + 0.1 * diff(range(c_df$Plasma, na.rm = TRUE)),
                color = Time_Point,
                label = paste0(recode(Time_Point, Pre = "Pre", T0 = "0h", T24 = "24h"),
                               " ρ = ", round(Correlation, 2)))) +
  labs(x = "Saliva MANSC1 (z-score)", y = "Plasma MANSC1 (z-score)") +
  theme_pubr(base_size = 6) + theme(text = element_text(size = 6), legend.position = "none") + theme_strokes
ggsave(here("figures", "figure_4", "c_mansc1_example.pdf"), pc, width = 8.2, height = 4.5, units = "cm", dpi = 600)
SHEETS[["c_MANSC1"]] <- c_df

# ---- panel f: per-protein cross-fluid Spearman rho, Pre vs 0 h ---------------
time_points <- c("Pre", "T0", "T24")
res_rows <- list()
for (core in matching_cores) {
  pr <- filtered_data[grep(paste0("^", core, "_"), protein_names), , drop = FALSE]
  if (nrow(pr) < 2) next
  for (tm in time_points) {
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
results <- bind_rows(res_rows) %>% mutate(Time_Point = recode(Time_Point, T0 = "0h", T24 = "24h"))
corr_data <- results %>%
  filter(!Source %in% c("urine_urine", "saliva_saliva")) %>%
  mutate(Source = str_replace(Source, "olink", "plasma")) %>%
  group_by(Protein, Source, Time_Point) %>%
  summarise(Correlation = mean(Correlation, na.rm = TRUE), .groups = "drop") %>%
  filter(Time_Point %in% c("Pre", "0h")) %>%
  pivot_wider(id_cols = c(Protein, Source), names_from = Time_Point,
              values_from = Correlation, names_prefix = "Time_") %>%
  drop_na(Time_Pre, Time_0h)
f_pal <- c("Both <-0.5" = "#1E88E5", "Both >0.5" = "#1E88E5", "Top-Left" = "#FFA07A",
           "Bottom-Right" = "#FFA07A", ">0.5 or <-0.5 and ~0" = "#004D40",
           "~0 and >0.5 or <-0.5" = "#D81B60", "Unlabeled" = "grey80")
corr_data <- corr_data %>% mutate(
  Label = case_when(
    Time_Pre < -0.5 & Time_0h < -0.5 ~ Protein, Time_Pre > 0.5 & Time_0h > 0.5 ~ Protein,
    Time_Pre > 0.4 & Time_0h < -0.4 ~ Protein, Time_Pre < -0.4 & Time_0h > 0.4 ~ Protein,
    (Time_Pre > 0.5 | Time_Pre < -0.5) & Time_0h > -0.1 & Time_0h < 0.1 ~ Protein,
    (Time_0h > 0.5 | Time_0h < -0.5) & Time_Pre > -0.1 & Time_Pre < 0.1 ~ Protein,
    TRUE ~ NA_character_),
  Color = case_when(
    Time_Pre < -0.5 & Time_0h < -0.5 ~ "Both <-0.5", Time_Pre > 0.5 & Time_0h > 0.5 ~ "Both >0.5",
    Time_Pre > 0.4 & Time_0h < -0.4 ~ "Top-Left", Time_Pre < -0.4 & Time_0h > 0.4 ~ "Bottom-Right",
    (Time_Pre > 0.5 | Time_Pre < -0.5) & Time_0h > -0.1 & Time_0h < 0.1 ~ ">0.5 or <-0.5 and ~0",
    (Time_0h > 0.5 | Time_0h < -0.5) & Time_Pre > -0.1 & Time_Pre < 0.1 ~ "~0 and >0.5 or <-0.5",
    TRUE ~ "Unlabeled"))
f_plots <- corr_data %>% split(.$Source) %>% lapply(function(df) {
  ggplot(df, aes(Time_Pre, Time_0h, color = Color)) +
    geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black", size = 0.5) +
    geom_point(size = 1, alpha = 0.7) +
    ggrepel::geom_text_repel(aes(label = Label), size = 2, force = 20, max.overlaps = Inf, na.rm = TRUE) +
    scale_color_manual(values = f_pal) +
    labs(title = paste("Across:", unique(df$Source)), x = expression(rho ~ "at Pre"), y = expression(rho ~ "at 0h")) +
    theme_pubr(base_size = 6) + theme(text = element_text(size = 6), legend.position = "none") + theme_strokes
})
pf <- wrap_plots(f_plots, ncol = 3)
ggsave(here("figures", "figure_4", "f_rho_pre_vs_0h.pdf"), pf, width = 16.5, height = 5.5, units = "cm", dpi = 600)
SHEETS[["f_rho_pre_vs_0h"]] <- as.data.frame(corr_data)

# =============================================================================
# panels d & e — cross-fluid networks at 0 h and 1 h (ORA-aware blobs)
# plotting lives in R/figures/helpers_crossfluid_network.R (shared with ED1).
# =============================================================================
g0 <- cf_plot_network("0", here("figures", "figure_4", "d_network_0h.pdf"), NET_DIR)
g1 <- cf_plot_network("1", here("figures", "figure_4", "e_network_1h.pdf"), NET_DIR)
SHEETS[["d_network_0h_inset"]] <- cf_top_hub_inset(g0, here("figures", "figure_4", "d_network_0h_inset.pdf"))
SHEETS[["e_network_1h_inset"]] <- cf_top_hub_inset(g1, here("figures", "figure_4", "e_network_1h_inset.pdf"))
SHEETS[["d_network_0h_edges"]] <- igraph::as_data_frame(g0, "edges")
SHEETS[["e_network_1h_edges"]] <- igraph::as_data_frame(g1, "edges")

# Source Data for this figure is built separately by
# R/source_data/figure4_source_data.R (writes source_data/SourceData_Figure4.xlsx).
cat("Figure 4 panels written:", paste(names(SHEETS), collapse = ", "), "\n")
