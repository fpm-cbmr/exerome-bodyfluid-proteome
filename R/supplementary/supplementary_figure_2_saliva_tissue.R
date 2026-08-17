# =============================================================================
# Supplementary Figure 2 — Tissue and cell-type-specific contribution to the exercise-regulated salivary proteome (related to Fig. 2).
#
#   a  saliva proteins binned into HPA tissues (shape = first up/down, size = |beta|,
#      colour = temporal cluster)
#   b  Fisher's exact tissue enrichment across the 601 regulated saliva proteins
#   c  Fisher's exact cell-type enrichment across the 601 regulated saliva proteins
#
# Input:  data/top_tissue_per_protein_saliva.rda, data/sig_vs_tissue_cnt_saliva.rda,
#         data/sig_vs_celltype_cnt_saliva.rda
# Output: figures/supplementary/S2_*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(ggplot2); library(ggpubr); library(ggrepel)
})
OUT <- here("figures", "supplementary"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cluster_color_palette <- c("0"="#B0B0B0","1"="#1F78B4","2"="#33A02C","3"="#E31A1C",
                           "4"="#FDBF6F","5"="#FF7F00","6"="#6A3D9A","7"="#CAB2D6")

# ---- a: tissue-binned scatter -----------------------------------------------
load(here("data", "top_tissue_per_protein_saliva.rda"))
top_tissue_per_protein <- top_tissue_per_protein_saliva %>%
  mutate(b_high_numeric = abs(as.numeric(b.high)))
tissue_order <- top_tissue_per_protein %>% group_by(tissue) %>% summarise(Count = n()) %>%
  arrange(desc(Count)) %>% .$tissue
top_tissue_per_protein$tissue <- factor(top_tissue_per_protein$tissue, levels = tissue_order)
top_tissue_per_protein$first_sig <- factor(top_tissue_per_protein$first_sig,
                                           levels = c("0", "0.5", "1", "3", "24", "NA"))
top_tissue_per_protein <- top_tissue_per_protein %>%
  filter(!grepl("^NA.", Assay), !is.na(b_high_dir))

top_assays <- top_tissue_per_protein %>% arrange(desc(-log10(fdr.aov))) %>%
  group_by(tissue) %>% slice_max(-log10(fdr.aov), n = 1) %>% ungroup()
max_value <- 16
pa <- ggplot(top_tissue_per_protein, aes(x = tissue, y = -log10(fdr.aov))) +
  geom_point(aes(color = factor(cluster), size = b_high_numeric, shape = b_high_dir), alpha = 0.8) +
  scale_shape_manual(values = c("up" = 24, "down" = 25)) +
  scale_size_continuous(range = c(0.5, 3), breaks = c(0.5, 1.5, 2.5)) +
  scale_color_manual(values = cluster_color_palette) +
  labs(title = NULL, x = NULL, y = expression(-log[10]("q-value"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "top", legend.direction = "horizontal", text = element_text(size = 6),
        legend.title = element_text(size = 6), legend.text = element_text(size = 6),
        legend.box = "horizontal", legend.key.size = unit(2, "mm"),
        panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = NA)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(data = top_assays, aes(label = Assay), size = 2, max.overlaps = Inf, force = 20,
                  segment.color = "grey50", segment.size = 0.5, direction = "y", nudge_y = 0.3) +
  ylim(c(0, max_value)) +
  guides(color = guide_legend(title = "cluster", nrow = 2, byrow = TRUE),
         size = guide_legend(title = "beta", nrow = 1, byrow = TRUE),
         shape = guide_legend(title = "first_regulation", nrow = 2, byrow = TRUE))
ggsave(file.path(OUT, "S2a_saliva_tissue_summary.pdf"), pa, width = 16.5, height = 7, units = "cm", dpi = 600)

# ---- Fisher dotplot for a tissue/cell-type count table ----------------------
enrichment_dotplot <- function(df, feat, title, low_col) {
  df <- df %>% arrange(padj)
  df$.feat <- df[[feat]]
  ordered <- reorder(df$.feat, df$fisher_p)
  last_sig <- levels(ordered)[max(which(df$padj < 0.05))]
  vline_x  <- which(levels(ordered) == last_sig) + 1
  ggplot(df, aes(x = reorder(.feat, fisher_p), y = fisher_OR, fill = -log10(fisher_p))) +
    geom_point(shape = 21, aes(size = sig + not_sig, stroke = ifelse(padj < 0.05, 1, 0))) +
    geom_vline(xintercept = vline_x, linetype = "dashed", color = "red") +
    scale_fill_gradient(high = "lightyellow", low = low_col, name = "-log10(p-value)",
                        trans = "reverse") +
    scale_size_continuous(name = "Proteins", range = c(0.5, 3)) +
    labs(title = title, x = NULL, y = "Odds ratio") +
    theme_bw() +
    theme(legend.position = c(0.6, 1), legend.direction = "horizontal",
          legend.background = element_blank(), legend.key = element_rect(fill = "transparent"),
          legend.key.size = unit(2, "mm"), legend.title = element_text(size = 6),
          legend.text = element_text(size = 6), legend.justification = "top",
          axis.title = element_text(size = 6),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6), text = element_text(size = 6),
          panel.grid = element_blank()) +
    guides(fill = guide_colorbar(order = 1), size = guide_legend(order = 2),
           stroke = guide_legend(order = 3, override.aes = list(size = 4)))
}

# ---- b: tissue enrichment ---------------------------------------------------
load(here("data", "sig_vs_tissue_cnt_saliva.rda"))
pb <- enrichment_dotplot(sig_vs_tissue_cnt_saliva, "tissue", "Saliva - Tissue Enrichment", "#FFA07A")
ggsave(file.path(OUT, "S2b_saliva_tissue_enrichment.pdf"), pb, width = 16.5, height = 6.5, units = "cm", dpi = 600)

# ---- c: cell-type enrichment ------------------------------------------------
load(here("data", "sig_vs_celltype_cnt_saliva.rda"))
pc <- enrichment_dotplot(sig_vs_celltype_cnt_saliva, "celltype", "Saliva - Cell Type Enrichment", "#FFA07A")
ggsave(file.path(OUT, "S2c_saliva_celltype_enrichment.pdf"), pc, width = 16.5, height = 6.5, units = "cm", dpi = 600)

# source data now built by R/source_data/ (per-figure script)
cat("Supplementary Figure 2 (saliva tissue/cell-type) written to figures/supplementary/S2_*.pdf\n")
