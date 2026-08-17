# =============================================================================
# Figure 3 — The exercise-induced plasma proteome reflects multi-tissue origins and disease associations.
#
# Panels:
#   a plasma t-SNE  b baseline-vs-residual variance scatter  c NADK/FOLR3 examples
#   d plasma cluster dynamics  e plasma cluster ORA  f secretome proportion
#   g tissue  h cell-type  i disease (HR<1)
#
# Inputs (produced by the analysis pipeline):
#   data/plasma_tsne.rda                        <- 04_plasma_tsne.R                  (a)
#   data/res.olink.linear.rda, data/olink.exerome.dat.rda,
#     data/plasma_variance_decomposition.rda    <- 04_plasma_lmm_clusters.R          (a-d,f,g,h)
#   data/top_tissue_per_protein_olink.rda,
#     data/cluster_vs_{tissue,cell}_cnt_olink.rda <- 08_plasma_hpa_enrichment.R       (f,g,h)
#   data/res.enrich.olink.cluster.rda           <- 07_plasma_cluster_ora.R          (e)
#   doc/supplemental_tables/disease_enrichment_decreasing_plasma.csv  (disease, HR<1) (i)
#
# Output: figures/figure_3/*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2); library(ggpubr)
  library(writexl); library(Polychrome); library(ComplexHeatmap); library(circlize); library(grid)
  library(readr); library(stringr); library(ggrepel)
})
source(here("R/package_loading.R")); source(here("R/functions_loading.R")); source(here("R/figure_defaults.R"))
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate; rename <- dplyr::rename
group_by <- dplyr::group_by; summarise <- dplyr::summarise; arrange <- dplyr::arrange
count <- dplyr::count; slice_head <- dplyr::slice_head
dir.create(here("figures", "figure_3"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()
cluster_color_palette <- c("0"="#B0B0B0","1"="#1F78B4","2"="#33A02C","3"="#E31A1C",
                           "4"="#FDBF6F","5"="#FF7F00","6"="#6A3D9A","7"="#CAB2D6")
participant_palette <- function(part) { lv <- sort(unique(as.character(part)))
  setNames(palette36.colors(max(25, length(lv)))[seq_along(lv)], lv) }

load(here("data", "res.olink.linear.rda"))
load(here("data", "olink.exerome.dat.rda"))

# ---- panel a: plasma t-SNE ---------------------------------------------------
load(here("data", "plasma_tsne.rda"))
pal_p <- participant_palette(plasma_tsne$participant)
plasma_tsne$participant <- factor(as.character(plasma_tsne$participant), levels = names(pal_p))
pa <- ggplot(plasma_tsne, aes(tSNE1, tSNE2, color = participant)) +
  geom_point(size = 1.5) + scale_color_manual(values = pal_p) +
  labs(title = "Plasma Proteome Individuality", x = "t-SNE1", y = "t-SNE2") +
  theme_bw(base_size = 6) +
  theme(legend.position = "none", panel.grid = element_blank(),
        text = element_text(size = 6), plot.title = element_text(size = 6),
        axis.title = element_text(size = 6), axis.text = element_text(size = 6))
ggsave(here("figures", "figure_3", "a_plasma_tsne.pdf"), pa, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["a_plasma_tsne"]] <- plasma_tsne

# ---- panel d: plasma temporal cluster dynamics -------------------------------
idcol <- "OlinkID"
merged_p <- olink.exerome.dat %>%
  pivot_longer(cols = -c(time_label, subject, replicate, gender, t.factor),
               names_to = idcol, values_to = "value") %>%
  left_join(res.olink.linear, by = idcol) %>% filter(!is.na(cluster))
summary_p <- merged_p %>% group_by(time_label, cluster) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = ifelse(n() > 1, sd(value, na.rm = TRUE), 0), n = n(), .groups = "drop") %>%
  mutate(se = ifelse(n > 1, sd_value / sqrt(n), 0), ci_lower = mean_value - 1.96 * se, ci_upper = mean_value + 1.96 * se)
nprot_p <- merged_p %>% group_by(cluster) %>% summarise(count = n_distinct(.data[[idcol]]), .groups = "drop") %>%
  mutate(label = paste0("C", cluster, ": ", count), color = cluster_color_palette[as.character(cluster)])
pd <- ggplot(summary_p, aes(time_label, mean_value, color = factor(cluster), group = factor(cluster))) +
  geom_rect(aes(xmin = 0, xmax = 4, ymin = -Inf, ymax = Inf), fill = "grey95", color = NA, alpha = 0.7) +
  geom_line() + geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = factor(cluster)), alpha = 0.2, color = NA) +
  scale_color_manual(values = cluster_color_palette) + scale_fill_manual(values = cluster_color_palette) +
  geom_hline(yintercept = 0, color = "black") + scale_x_continuous(breaks = c(-1, 0, 1, 3, 24)) +
  labs(x = "Time [h]", y = "Plasma Mean Z-score") + theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(), legend.position = c(1, 1.03), legend.justification = c("right", "top"),
        legend.background = element_blank(), legend.key.size = unit(2, "mm")) + theme_strokes
ylp <- range(-1, summary_p$ci_upper)
yposp <- seq(ylp[1] + 0.15 * diff(ylp), ylp[1] - 0.2 * diff(ylp), length.out = nrow(nprot_p))
for (i in seq_len(nrow(nprot_p)))
  pd <- pd + annotate("text", x = max(summary_p$time_label) - 3.5, y = yposp[i],
                      label = nprot_p$label[i], color = nprot_p$color[i], hjust = 0, size = 2)
ggsave(here("figures", "figure_3", "d_plasma_dynamics.pdf"), pd, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["d_plasma_dynamics"]] <- summary_p

# ---- panels g,h: tissue & cell-type enrichment heatmaps ----------------------
enrich_heatmap <- function(cnt, row_key, out_pdf, title) {
  cnt <- cnt %>% filter(!is.na(.data[[row_key]]), .data[[row_key]] != "") %>%
    distinct(.data[[row_key]], cluster, .keep_all = TRUE)
  keep_rows <- cnt %>% group_by(.data[[row_key]]) %>%
    summarise(any_sig = any(padj < 0.05, na.rm = TRUE), .groups = "drop") %>% filter(any_sig) %>% pull(1)
  cnt2 <- cnt %>% filter(.data[[row_key]] %in% keep_rows)   # keep cluster 0 (not-regulated column)
  or <- cnt2 %>% dplyr::select(all_of(row_key), cluster, fisher_OR) %>%
    tidyr::pivot_wider(names_from = cluster, values_from = fisher_OR) %>% as.data.frame()
  rownames(or) <- or[[row_key]]; or[[row_key]] <- NULL
  pmat <- cnt2 %>% dplyr::select(all_of(row_key), cluster, padj) %>%
    tidyr::pivot_wider(names_from = cluster, values_from = padj) %>% as.data.frame()
  rownames(pmat) <- pmat[[row_key]]; pmat[[row_key]] <- NULL
  m <- as.matrix(or); m[m == 0] <- NA; m <- log2(m); pm <- as.matrix(pmat)[rownames(m), colnames(m)]
  colnames(m) <- paste0("C", colnames(m))
  m_cl <- m; m_cl[!is.finite(m_cl)] <- 0   # NA-safe copy for row clustering (cluster 0 has many empty tissues)
  ht <- Heatmap(m, name = "Log2(OR)", col = colorRamp2(c(-2, 0, 2), c("#006ae3", "white", "#d70000")),
                na_col = "grey", cluster_columns = FALSE,
                cluster_rows = as.dendrogram(hclust(dist(m_cl))), column_title = title,
                column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize = 6), column_title_gp = gpar(fontsize = 6),
                cell_fun = function(j, i, x, y, w, h, fill) if (!is.na(pm[i, j]) && pm[i, j] < 0.05) grid.text("*", x, y, gp = gpar(fontsize = 8)),
                heatmap_legend_param = list(labels_gp = gpar(fontsize = 6), title_gp = gpar(fontsize = 6)))
  pdf(out_pdf, width = 5.5 / 2.54, height = 5 / 2.54); draw(ht); dev.off()
  or2 <- or; or2[[row_key]] <- rownames(or); or2
}
load(here("data", "cluster_vs_tissue_cnt_olink.rda"))
SHEETS[["g_plasma_tissue"]] <- enrich_heatmap(cluster_vs_tissue_cnt_olink, "tissue",
  here("figures", "figure_3", "g_plasma_tissue.pdf"), "Plasma - Tissue")
load(here("data", "cluster_vs_cell_cnt_olink.rda"))
SHEETS[["h_plasma_celltype"]] <- enrich_heatmap(cluster_vs_cell_cnt_olink, "celltype",
  here("figures", "figure_3", "h_plasma_celltype.pdf"), "Plasma - Cell type")

# ---- panel f: plasma secretome proportion ------------------------------------
load(here("data", "top_tissue_per_protein_olink.rda"))
tt <- top_tissue_per_protein_olink %>% dplyr::select(OlinkID, Secretome.function) %>%
  left_join(res.olink.linear %>% dplyr::select(OlinkID, cluster), by = "OlinkID")  # correct clusters
prop_p <- tt %>% mutate(cluster = tidyr::replace_na(cluster, 0)) %>%   # NA cluster = not regulated (cluster 0)
  mutate(Secretome.function = na_if(Secretome.function, "")) %>%
  filter(!is.na(Secretome.function)) %>%
  group_by(Secretome.function, cluster) %>% summarise(count = n(), .groups = "drop") %>%
  group_by(Secretome.function) %>% mutate(total = sum(count), proportion = count / total) %>% ungroup()
sec_ord_p <- prop_p %>% group_by(Secretome.function) %>% summarise(tc = sum(count), .groups = "drop") %>% arrange(desc(tc)) %>% pull(Secretome.function)
prop_p$Secretome.function <- factor(prop_p$Secretome.function, levels = sec_ord_p)
pf <- ggplot(prop_p, aes(Secretome.function, proportion, fill = factor(cluster))) +
  geom_bar(stat = "identity", position = "stack") + labs(x = "Secretome Function", y = "Proportion") +
  theme_bw(base_size = 6) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none",
        panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = NA)) +
  scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values = cluster_color_palette) + theme_strokes
ggsave(here("figures", "figure_3", "f_plasma_secretome.pdf"), pf, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["f_plasma_secretome"]] <- prop_p

# ---- panel i: plasma cluster disease enrichment (HR<1) -----------------------
res_dec <- read_csv(here("doc", "supplemental_tables", "disease_enrichment_decreasing_plasma.csv"),
                    show_col_types = FALSE) %>%
  mutate(
    Cluster = as.character(Cluster),
    Odds_Ratio = as.numeric(Odds_Ratio),
    P_value = as.numeric(P_value),
    padj = as.numeric(padj)
  )

if (!"Log10P_signed" %in% names(res_dec)) {
  res_dec <- res_dec %>% mutate(Log10P_signed = if_else(Odds_Ratio > 1, -log10(P_value), log10(P_value)))
} else {
  res_dec <- res_dec %>% mutate(Log10P_signed = as.numeric(Log10P_signed))
}

if (!"Point_Color" %in% names(res_dec)) {
  res_dec <- res_dec %>% mutate(Point_Color = if_else(padj < 0.05, Cluster, "Non-significant"))
}

res_dec <- res_dec %>%
  mutate(Point_Color = as.character(Point_Color)) %>%
  filter(Cluster != "0")
disease_order_p <- res_dec %>% group_by(Disease) %>%
  summarise(m = max(Log10P_signed, na.rm = TRUE), .groups = "drop") %>% arrange(desc(m)) %>% pull(Disease)
res_dec$Disease <- factor(res_dec$Disease, levels = disease_order_p)   # sort by -log10 P
top5_p <- res_dec %>% filter(padj < 0.05) %>% group_by(Cluster) %>% arrange(P_value) %>% slice_head(n = 5) %>% ungroup()
pi <- ggplot(res_dec, aes(Disease, Log10P_signed)) +
  geom_point(aes(color = Point_Color, size = Odds_Ratio), alpha = 0.7) +
  scale_color_manual(values = c(cluster_color_palette, "Non-significant" = "grey90"), name = "Cluster") +
  scale_size_continuous(range = c(0.5, 3), name = "Odds Ratio") +
  geom_text_repel(data = top5_p, aes(Disease, Log10P_signed, label = Disease, color = Cluster),
                  inherit.aes = FALSE, size = 2, force = 40, max.overlaps = 50, show.legend = FALSE) +
  labs(title = "Disease Enrichment Across Plasma Clusters (HR<1)", x = "Disease", y = "Signed -log10(P-value)") +
  theme_pubr(base_size = 6) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        text = element_text(size = 6), legend.position = "right", legend.key.size = unit(2, "mm")) + theme_strokes
ggsave(here("figures", "figure_3", "i_plasma_disease.pdf"), pi, width = 11, height = 6, units = "cm", dpi = 600)
SHEETS[["i_plasma_disease"]] <- res_dec %>%
  mutate(Odds_Ratio = as.character(Odds_Ratio), Cluster = as.character(Cluster))

# ---- panel b: baseline vs residual variance scatter --------------------------
# variance decomposition from R/analysis/04_plasma_lmm_clusters.R (not random forest)
load(here("data", "plasma_variance_decomposition.rda"))
vdf <- plasma_variance_decomposition %>% dplyr::select(OlinkID, Assay, Baseline_Variance, Residual_Variance) %>%
  left_join(res.olink.linear %>% dplyr::select(OlinkID, cluster), by = "OlinkID")
vc0 <- vdf %>% filter(is.na(cluster) | cluster == 0)
vco <- vdf %>% filter(!is.na(cluster) & cluster != 0)
pb <- ggplot() +
  geom_point(data = vc0, aes(Residual_Variance, Baseline_Variance), color = "grey85", size = 0.5) +
  geom_point(data = vco, aes(Residual_Variance, Baseline_Variance, color = factor(cluster)), size = 0.5) +
  geom_text_repel(data = vdf %>% filter(Assay %in% c("NADK", "FOLR3")),
                  aes(Residual_Variance, Baseline_Variance, label = Assay), size = 2, max.overlaps = Inf) +
  scale_color_manual(values = cluster_color_palette) +
  labs(x = "Residual variance", y = "Baseline variance") +
  theme_pubr(base_size = 6) + theme(legend.position = "none", text = element_text(size = 6)) + theme_strokes
ggsave(here("figures", "figure_3", "b_variance_scatter.pdf"), pb, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["b_variance_scatter"]] <- vdf

# ---- panel c: NPX over time for NADK and FOLR3 -------------------------------
ex_ids <- res.olink.linear %>% filter(Assay %in% c("NADK", "FOLR3")) %>% dplyr::select(OlinkID, Assay)
pc_df <- olink.exerome.dat %>% dplyr::select(time_label, replicate, all_of(ex_ids$OlinkID)) %>%
  pivot_longer(all_of(ex_ids$OlinkID), names_to = "OlinkID", values_to = "NPX") %>%
  left_join(ex_ids, by = "OlinkID")
pcp <- ggplot(pc_df, aes(time_label, NPX, color = factor(replicate), group = interaction(replicate, Assay))) +
  geom_line(alpha = 0.6) + facet_wrap(~ Assay, scales = "free_y") +
  scale_color_manual(values = participant_palette(pc_df$replicate)) +
  scale_x_continuous(breaks = c(-1, 0, 1, 3, 24)) + labs(x = "Time [h]", y = "NPX") +
  theme_pubr(base_size = 6) + theme(legend.position = "none", text = element_text(size = 6),
        strip.text = element_text(size = 6)) + theme_strokes
ggsave(here("figures", "figure_3", "c_individual_proteins.pdf"), pcp, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["c_individual_proteins"]] <- pc_df

# ---- panel e: plasma selected cluster ORA (GO:BP) ----------------------------
load(here("data", "res.enrich.olink.cluster.rda"))
specific_terms_p <- c("Fat digestion and absorption", "leukocyte migration", "response to hormone",
  "neuron projection development", "neurotransmitter secretion", "response to auditory stimulus",
  "neutrophil degranulation", "Degradation of the extracellular matrix", "tissue remodeling",
  "vesicle-mediated transport", "Signaling by Rho GTPases", "endosomal transport",
  "regulation of proteolysis", "apoptotic signaling pathway", "Cytokine Signaling in Immune system",
  "angiogenesis", "regulation of transmembrane transport", "Muscle contraction")
wrap_text <- function(x, width = 35) sapply(x, function(y) paste(strwrap(y, width = width), collapse = "\n"))
plot_data_p <- res.enrich.olink.cluster %>%
  filter(!stringr::str_detect(result.term_name, "root term")) %>%
  group_by(cluster, result.intersection_size) %>% dplyr::slice_min(result.p_value, n = 1, with_ties = FALSE) %>% ungroup() %>%
  filter(result.term_name %in% specific_terms_p) %>% mutate(cluster.paper = paste0("C", cluster))
plot_data_p$result.term_name <- wrap_text(plot_data_p$result.term_name)
pe <- ggplot(plot_data_p, aes(reorder(result.term_name, -result.p_value), -log10(result.p_value), fill = as.factor(cluster))) +
  geom_bar(stat = "identity", width = 0.5) + labs(x = NULL, y = "-log10(adj.p-val)") +
  scale_fill_manual(values = cluster_color_palette) + theme_pubr(base_size = 6) +
  theme(axis.text.x = element_text(size = 6, angle = 70, hjust = 1, vjust = 1), legend.position = "none", axis.ticks = element_blank()) +
  facet_wrap(~ cluster.paper, scales = "free_x", nrow = 1) + theme_strokes
ggsave(here("figures", "figure_3", "e_plasma_cluster_ORA.pdf"), pe, width = 11, height = 5, units = "cm", dpi = 600)
SHEETS[["e_plasma_cluster_ORA"]] <- plot_data_p

# Source Data for this figure is built separately by
# R/source_data/figure3_source_data.R (writes source_data/SourceData_Figure3.xlsx).
cat("Figure 3 panels:", paste(names(SHEETS), collapse = ", "), "\n")
