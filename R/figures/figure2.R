# =============================================================================
# Figure 2 — Saliva and urine proteomes capture compartment-specific responses.
#
# Panels:
#   a saliva t-SNE  b saliva cluster dynamics  c saliva secretome proportion
#   d saliva cluster ORA  e saliva tissue  f saliva cell-type  g saliva disease (HR>1)
#   h urine t-SNE  i urine cluster dynamics  j urine cell-type
#
# Inputs (produced by the analysis pipeline):
#   data/saliva_tsne.rda / data/urine_tsne.rda            <- 02_saliva_tsne.R / 03_urine_tsne.R   (a,h)
#   data/res.saliva.linear.rda, data/exerome.dat.saliva.rda <- 02_saliva_lmm_clusters.R           (b)
#   data/res.urine.linear.rda,  data/exerome.dat.urine.rda  <- 03_urine_lmm_clusters.R            (i)
#   data/top_tissue_per_protein_saliva.rda,
#     data/cluster_vs_{tissue,cell}_cnt_saliva.rda         <- 05_saliva_hpa_enrichment.R          (c,e,f)
#   data/cluster_vs_cell_cnt_urine.rda                     <- 09_urine_hpa_enrichment.R           (j)
#   data/res.enrich.saliva.cluster.rda                     <- 06_saliva_cluster_ora.R            (d)
#   doc/supplemental_tables/disease_enrichment_increasing.csv  (saliva cluster disease enrichment)  (g)
#
# Output: figures/figure_2/*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2); library(ggpubr)
  library(writexl); library(Polychrome)
})
# participant color coding used throughout: 25 distinct Polychrome colors, sorted
participant_palette <- function(part) {
  lv <- sort(unique(as.character(part)))
  setNames(palette36.colors(max(25, length(lv)))[seq_along(lv)], lv)
}
source(here("R/package_loading.R"))
source(here("R/functions_loading.R"))
source(here("R/figure_defaults.R"))           # 0.5 pt strokes
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
rename <- dplyr::rename; group_by <- dplyr::group_by; summarise <- dplyr::summarise
arrange <- dplyr::arrange; count <- dplyr::count; slice_head <- dplyr::slice_head
dir.create(here("figures", "figure_2"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("source_data"),         showWarnings = FALSE)
SHEETS <- list()   # accumulate source-data sheets

# ---- panel a: saliva t-SNE (seeded; from analysis/02_saliva_tsne.R) ----------
load(here("data", "saliva_tsne.rda"))
pal_sal <- participant_palette(saliva_tsne$participant)
saliva_tsne$participant <- factor(as.character(saliva_tsne$participant), levels = names(pal_sal))
pa <- ggplot(saliva_tsne, aes(tSNE1, tSNE2, color = participant)) +
  geom_point(size = 1.5) + scale_color_manual(values = pal_sal) +
  labs(title = "Saliva Proteome Individuality", x = "t-SNE1", y = "t-SNE2") +
  theme_bw(base_size = 6) +
  theme(legend.position = "none", panel.grid = element_blank(),
        text = element_text(size = 6), plot.title = element_text(size = 6),
        axis.title = element_text(size = 6), axis.text = element_text(size = 6))
ggsave(here("figures", "figure_2", "a_saliva_tsne.pdf"), pa, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["a_saliva_tsne"]] <- saliva_tsne

# ---- panel h: urine t-SNE (seeded; from analysis/03_urine_tsne.R) ------------
load(here("data", "urine_tsne.rda"))
pal_uri <- participant_palette(urine_tsne$participant)
urine_tsne$participant <- factor(as.character(urine_tsne$participant), levels = names(pal_uri))
ph <- ggplot(urine_tsne, aes(tSNE1, tSNE2, color = participant)) +
  geom_point(size = 1.5) + scale_color_manual(values = pal_uri) +
  labs(title = "Urine Proteome Individuality", x = "t-SNE1", y = "t-SNE2") +
  theme_bw(base_size = 6) +
  theme(legend.position = "none", panel.grid = element_blank(),
        text = element_text(size = 6), plot.title = element_text(size = 6),
        axis.title = element_text(size = 6), axis.text = element_text(size = 6))
ggsave(here("figures", "figure_2", "h_urine_tsne.pdf"), ph, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["h_urine_tsne"]] <- urine_tsne

# ---- panel b: saliva temporal cluster dynamics (z-scores by cluster) ---------
cluster_color_palette <- c("0"="#B0B0B0","1"="#1F78B4","2"="#33A02C","3"="#E31A1C",
                           "4"="#FDBF6F","5"="#FF7F00","6"="#6A3D9A","7"="#CAB2D6")
load(here("data", "res.saliva.linear.rda"))
load(here("data", "exerome.dat.saliva.rda"))
merged_sal <- exerome.dat.saliva %>%
  pivot_longer(cols = -c(time_label, subject, replicate, gender, t.factor),
               names_to = "Assay", values_to = "value") %>%
  left_join(res.saliva.linear, by = "Assay") %>% filter(!is.na(cluster))
summary_sal <- merged_sal %>% group_by(time_label, cluster) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = ifelse(n() > 1, sd(value, na.rm = TRUE), 0), n = n(), .groups = "drop") %>%
  mutate(se = ifelse(n > 1, sd_value / sqrt(n), 0),
         ci_lower = mean_value - 1.96 * se, ci_upper = mean_value + 1.96 * se)
nprot_sal <- merged_sal %>% group_by(cluster) %>% summarise(count = n_distinct(Assay), .groups = "drop") %>%
  mutate(label = paste0("C", cluster, ": ", count), color = cluster_color_palette[as.character(cluster)])
pb <- ggplot(summary_sal, aes(time_label, mean_value, color = factor(cluster), group = factor(cluster))) +
  geom_rect(aes(xmin = 0, xmax = 4, ymin = -Inf, ymax = Inf), fill = "grey95", color = NA, alpha = 0.7) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = factor(cluster)), alpha = 0.2, color = NA) +
  scale_color_manual(values = cluster_color_palette) + scale_fill_manual(values = cluster_color_palette) +
  geom_hline(yintercept = 0, color = "black") +
  scale_x_continuous(breaks = c(-1, 0, 1, 3, 24)) +
  labs(x = "Time [h]", y = "Saliva Mean Z-score") +
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(), legend.position = c(1, 1.03),
        legend.justification = c("right", "top"), legend.background = element_blank(),
        legend.key.size = unit(2, "mm")) + theme_strokes
yl <- range(-1, summary_sal$ci_upper)
ypos <- seq(yl[1] + 0.15 * diff(yl), yl[1] - 0.2 * diff(yl), length.out = nrow(nprot_sal))
for (i in seq_len(nrow(nprot_sal)))
  pb <- pb + annotate("text", x = max(summary_sal$time_label) - 3.5, y = ypos[i],
                      label = nprot_sal$label[i], color = nprot_sal$color[i], hjust = 0, size = 2)
ggsave(here("figures", "figure_2", "b_saliva_dynamics.pdf"), pb, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["b_saliva_dynamics"]] <- summary_sal

# ---- panel c: saliva secreted-protein proportion by cluster ------------------
load(here("data", "top_tissue_per_protein_saliva.rda"))
prop_sal <- top_tissue_per_protein_saliva %>%
  mutate(cluster = tidyr::replace_na(cluster, 0)) %>%   # NA cluster = not regulated (cluster 0)
  mutate(Secretome.function = na_if(Secretome.function, "")) %>% filter(!is.na(Secretome.function)) %>%
  group_by(Secretome.function, cluster) %>% summarise(count = n(), .groups = "drop") %>%
  group_by(Secretome.function) %>% mutate(total = sum(count), proportion = count / total) %>% ungroup()
sec_order <- prop_sal %>% group_by(Secretome.function) %>% summarise(tc = sum(count), .groups = "drop") %>%
  arrange(desc(tc)) %>% pull(Secretome.function)
prop_sal$Secretome.function <- factor(prop_sal$Secretome.function, levels = sec_order)
pc <- ggplot(prop_sal, aes(Secretome.function, proportion, fill = factor(cluster))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Secretome Function", y = "Proportion", fill = "Cluster") +
  theme_bw(base_size = 6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none",
        panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = NA)) +
  scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values = cluster_color_palette) + theme_strokes
ggsave(here("figures", "figure_2", "c_saliva_secretome.pdf"), pc, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["c_saliva_secretome"]] <- prop_sal

# ---- panel d: saliva selected cluster ORA (GO:BP) ----------------------------
load(here("data", "res.enrich.saliva.cluster.rda"))
specific_terms <- list(
  "1" = c("sensory perception of taste", "skeletal system development", "peptide hormone processing", "collagen metabolic process"),
  "6" = c("response to external biotic stimulus", "regulation of vesicle-mediated transport", "programmed cell death"),
  "5" = c("aminoglycan metabolic process", "negative regulation of developmental growth", "Post-translational protein phosphorylation (REAC)"),
  "3" = c("Golgi-to-ER retrograde transport", "Autophagy"),
  "2" = c("Fatty acid degradation", "Glycolysis / Gluconeogenesis", "Lysine degradation"),
  "4" = c("Formation of the cornified envelope", "Keratinization"))
wrap_text <- function(x, width = 30) sapply(x, function(y) paste(strwrap(y, width = width), collapse = "\n"))
top_terms_sal <- res.enrich.saliva.cluster %>%
  filter(!stringr::str_detect(result.term_name, "root term")) %>%
  filter(paste0(cluster) %in% names(specific_terms)) %>%
  group_by(cluster) %>% filter(result.term_name %in% specific_terms[[as.character(unique(cluster))]]) %>% ungroup()
top_terms_sal$result.term_name <- wrap_text(top_terms_sal$result.term_name)
pd <- ggplot(top_terms_sal, aes(reorder(result.term_name, -result.p_value), -log10(result.p_value), fill = as.factor(cluster))) +
  geom_bar(stat = "identity", width = 0.5) + labs(x = NULL, y = "-log10(adj.p-val)") +
  scale_fill_manual(values = cluster_color_palette) + theme_pubr(base_size = 6) +
  theme(axis.text.x = element_text(size = 6, angle = 70, hjust = 1, vjust = 1), legend.position = "none",
        axis.ticks = element_blank()) +
  facet_wrap(~ cluster, scales = "free_x", nrow = 1, labeller = labeller(cluster = function(v) paste0("C", v))) +
  theme_strokes
ggsave(here("figures", "figure_2", "d_saliva_cluster_ORA.pdf"), pd, width = 11, height = 5, units = "cm", dpi = 600)
SHEETS[["d_saliva_cluster_ORA"]] <- top_terms_sal

# ---- panels e,f: tissue & cell-type enrichment heatmaps (Fisher log2 OR) ------
suppressMessages({library(ComplexHeatmap); library(circlize); library(grid)})
enrich_heatmap <- function(cnt, row_key, out_pdf, title) {
  cnt <- cnt %>% filter(!is.na(.data[[row_key]]), .data[[row_key]] != "") %>%
    distinct(.data[[row_key]], cluster, .keep_all = TRUE)
  keep_rows <- cnt %>% group_by(.data[[row_key]]) %>%
    summarise(any_sig = any(padj < 0.05, na.rm = TRUE), .groups = "drop") %>%
    filter(any_sig) %>% pull(1)
  cnt2 <- cnt %>% filter(.data[[row_key]] %in% keep_rows)   # keep cluster 0 (not-regulated column)
  or  <- cnt2 %>% dplyr::select(all_of(row_key), cluster, fisher_OR) %>%
    tidyr::pivot_wider(names_from = cluster, values_from = fisher_OR) %>% as.data.frame()
  rownames(or) <- or[[row_key]]; or[[row_key]] <- NULL
  pmat <- cnt2 %>% dplyr::select(all_of(row_key), cluster, padj) %>%
    tidyr::pivot_wider(names_from = cluster, values_from = padj) %>% as.data.frame()
  rownames(pmat) <- pmat[[row_key]]; pmat[[row_key]] <- NULL
  m <- as.matrix(or); m[m == 0] <- NA; m <- log2(m)
  pm <- as.matrix(pmat)[rownames(m), colnames(m)]
  colnames(m) <- paste0("C", colnames(m))
  m_cl <- m; m_cl[!is.finite(m_cl)] <- 0   # NA-safe copy for row clustering (cluster 0 has many empty tissues)
  ht <- Heatmap(m, name = "Log2(OR)", col = colorRamp2(c(-2, 0, 2), c("#006ae3", "white", "#d70000")),
                na_col = "grey", cluster_columns = FALSE,
                cluster_rows = as.dendrogram(hclust(dist(m_cl))), column_title = title,
                column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize = 6),
                column_title_gp = gpar(fontsize = 6),
                cell_fun = function(j, i, x, y, w, h, fill)
                  if (!is.na(pm[i, j]) && pm[i, j] < 0.05) grid.text("*", x, y, gp = gpar(fontsize = 8)),
                heatmap_legend_param = list(labels_gp = gpar(fontsize = 6), title_gp = gpar(fontsize = 6)))
  pdf(out_pdf, width = 5.5 / 2.54, height = 5 / 2.54); draw(ht); dev.off()
  or2 <- or; or2[[row_key]] <- rownames(or); or2
}
load(here("data", "cluster_vs_tissue_cnt_saliva.rda"))
SHEETS[["e_saliva_tissue"]] <- enrich_heatmap(cluster_vs_tissue_cnt_saliva, "tissue",
  here("figures", "figure_2", "e_saliva_tissue.pdf"), "Saliva - Tissue")
load(here("data", "cluster_vs_cell_cnt_saliva.rda"))
SHEETS[["f_saliva_celltype"]] <- enrich_heatmap(cluster_vs_cell_cnt_saliva, "celltype",
  here("figures", "figure_2", "f_saliva_celltype.pdf"), "Saliva - Cell type")


# ---- panel i: urine temporal cluster dynamics --------------------------------
load(here("data", "res.urine.linear.rda"))
load(here("data", "exerome.dat.urine.rda"))
merged_uri <- exerome.dat.urine %>%
  pivot_longer(cols = -c(time_label, subject, replicate, gender, t.factor),
               names_to = "Assay", values_to = "value") %>%
  left_join(res.urine.linear, by = "Assay") %>% filter(!is.na(cluster))
summary_uri <- merged_uri %>% group_by(time_label, cluster) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = ifelse(n() > 1, sd(value, na.rm = TRUE), 0), n = n(), .groups = "drop") %>%
  mutate(se = ifelse(n > 1, sd_value / sqrt(n), 0),
         ci_lower = mean_value - 1.96 * se, ci_upper = mean_value + 1.96 * se)
nprot_uri <- merged_uri %>% group_by(cluster) %>% summarise(count = n_distinct(Assay), .groups = "drop") %>%
  mutate(label = paste0("C", cluster, ": ", count), color = cluster_color_palette[as.character(cluster)])
pi <- ggplot(summary_uri, aes(time_label, mean_value, color = factor(cluster), group = factor(cluster))) +
  geom_rect(aes(xmin = 0, xmax = 4, ymin = -Inf, ymax = Inf), fill = "grey95", color = NA, alpha = 0.7) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = factor(cluster)), alpha = 0.2, color = NA) +
  scale_color_manual(values = cluster_color_palette) + scale_fill_manual(values = cluster_color_palette) +
  geom_hline(yintercept = 0, color = "black") + scale_x_continuous(breaks = c(-1, 0, 1, 3, 24)) +
  labs(x = "Time [h]", y = "Urine Mean Z-score") +
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(), legend.position = c(1, 1.03),
        legend.justification = c("right", "top"), legend.background = element_blank(),
        legend.key.size = unit(2, "mm")) + theme_strokes
yl2 <- range(-1, summary_uri$ci_upper)
ypos2 <- seq(yl2[1] + 0.15 * diff(yl2), yl2[1] - 0.2 * diff(yl2), length.out = nrow(nprot_uri))
for (i in seq_len(nrow(nprot_uri)))
  pi <- pi + annotate("text", x = max(summary_uri$time_label) - 3.5, y = ypos2[i],
                      label = nprot_uri$label[i], color = nprot_uri$color[i], hjust = 0, size = 2)
ggsave(here("figures", "figure_2", "i_urine_dynamics.pdf"), pi, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["i_urine_dynamics"]] <- summary_uri

# ---- panel j: urine cell-type enrichment heatmap -----------------------------
load(here("data", "cluster_vs_cell_cnt_urine.rda"))
SHEETS[["j_urine_celltype"]] <- enrich_heatmap(cluster_vs_cell_cnt_urine, "celltype",
  here("figures", "figure_2", "j_urine_celltype.pdf"), "Urine - Cell type")

# ---- panel g: saliva cluster disease enrichment (HR>1) -----------------------
suppressMessages({library(readr); library(ggrepel)})
res_inc <- read_csv(here("doc", "supplemental_tables", "disease_enrichment_increasing.csv"),
                    show_col_types = FALSE) %>%
  mutate(
    Cluster = as.character(Cluster),
    Odds_Ratio = as.numeric(Odds_Ratio),
    P_value = as.numeric(P_value),
    padj = as.numeric(padj),
    Log10P_signed = if_else(Odds_Ratio > 1, -log10(P_value), log10(P_value)),
    Point_Color = if_else(padj < 0.05, Cluster, "Non-significant")
  )
disease_order_sal <- res_inc %>% group_by(Disease) %>%
  summarise(m = max(Log10P_signed, na.rm = TRUE), .groups = "drop") %>% arrange(desc(m)) %>% pull(Disease)
res_inc$Disease <- factor(res_inc$Disease, levels = disease_order_sal)   # sort by -log10 P
top5_sal <- res_inc %>% filter(padj < 0.05) %>% group_by(Cluster) %>% arrange(P_value) %>% slice_head(n = 5) %>% ungroup()
pg <- ggplot(res_inc, aes(Disease, Log10P_signed)) +
  geom_point(aes(color = Point_Color, size = Odds_Ratio), alpha = 0.7) +
  scale_color_manual(values = c(cluster_color_palette, "Non-significant" = "grey90"), name = "Cluster") +
  scale_size_continuous(range = c(0.5, 3), name = "Odds Ratio") +
  geom_text_repel(data = top5_sal, aes(Disease, Log10P_signed, label = Disease, color = Cluster),
                  inherit.aes = FALSE, size = 2, force = 40, box.padding = 0.3, max.overlaps = 50, show.legend = FALSE) +
  labs(title = "Disease Enrichment Across Saliva Clusters (HR>1)", x = "Disease", y = "Signed -log10(P-value)") +
  theme_pubr(base_size = 6) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 6),
        legend.position = "right", legend.key.size = unit(2, "mm")) + theme_strokes
ggsave(here("figures", "figure_2", "g_saliva_disease.pdf"), pg, width = 11, height = 6, units = "cm", dpi = 600)
SHEETS[["g_saliva_disease"]] <- res_inc %>%
  mutate(Odds_Ratio = as.character(Odds_Ratio), Cluster = as.character(Cluster))

# Source Data for this figure is built separately by
# R/source_data/figure2_source_data.R (writes source_data/SourceData_Figure2.xlsx).
cat("Figure 2: panels ->", paste(names(SHEETS), collapse = ", "), "\n")
