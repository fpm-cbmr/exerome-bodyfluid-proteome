# =============================================================================
# Supplementary Figure 9 — Tissue and cell type-specific contribution to the exercise-regulated plasma proteome (related to Fig. 3).
#
#   a  plasma proteins binned into HPA tissues (shape = first up/down, size = |beta|,
#      colour = temporal cluster)
#   b  Fisher's exact tissue enrichment across the 2,156 regulated plasma proteins
#   c  Fisher's exact cell-type enrichment across the 2,156 regulated plasma proteins
#
# Input:  data/top_tissue_per_protein_olink.rda, data/res.olink.linear.rda,
#         data/sig_vs_tissue_cnt_olink.rda, data/sig_vs_celltype_cnt_olink.rda
# Output: figures/supplementary/S9_*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(ggplot2); library(ggpubr); library(ggrepel)
})
OUT <- here("figures", "supplementary"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cluster_color_palette <- c("0"="#B0B0B0","1"="#1F78B4","2"="#33A02C","3"="#E31A1C",
                           "4"="#FDBF6F","5"="#FF7F00","6"="#6A3D9A")

# ---- a: tissue-binned scatter (join fdr.aov / cluster / peak |beta| from LMM) --
load(here("data", "res.olink.linear.rda"))
beta_cols <- grep("^beta\\.exposure", names(res.olink.linear), value = TRUE)
lmm <- res.olink.linear %>%
  mutate(b_high_numeric = do.call(pmax, c(lapply(across(all_of(beta_cols)), abs), na.rm = TRUE))) %>%
  select(Assay, cluster, fdr.aov, b_high_numeric) %>% distinct(Assay, .keep_all = TRUE)

load(here("data", "top_tissue_per_protein_olink.rda"))
top_tissue_per_protein <- top_tissue_per_protein_olink %>%
  left_join(lmm, by = "Assay") %>%
  mutate(cluster = ifelse(is.na(cluster), 0, cluster)) %>%
  filter(!grepl("^NA.", Assay), !is.na(b_high_dir), !is.na(fdr.aov))
tissue_order <- top_tissue_per_protein %>% group_by(tissue) %>% summarise(Count = n()) %>%
  arrange(desc(Count)) %>% .$tissue
top_tissue_per_protein$tissue <- factor(top_tissue_per_protein$tissue, levels = tissue_order)
top_assays <- top_tissue_per_protein %>% arrange(desc(-log10(fdr.aov))) %>%
  group_by(tissue) %>% slice_max(-log10(fdr.aov), n = 1) %>% ungroup()
pa <- ggplot(top_tissue_per_protein, aes(x = tissue, y = -log10(fdr.aov))) +
  geom_point(aes(color = factor(cluster), size = b_high_numeric, shape = b_high_dir), alpha = 0.8) +
  scale_shape_manual(values = c("up" = 24, "down" = 25)) +
  scale_size_continuous(range = c(0.5, 3)) +
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
  guides(color = guide_legend(title = "cluster", nrow = 2, byrow = TRUE),
         size = guide_legend(title = "beta", nrow = 1, byrow = TRUE),
         shape = guide_legend(title = "first_regulation", nrow = 2, byrow = TRUE))
ggsave(file.path(OUT, "S9a_plasma_tissue_summary.pdf"), pa, width = 16.5, height = 7, units = "cm", dpi = 600)

# ---- Fisher dotplot helper (shared with S2/S4) ------------------------------
enrichment_dotplot <- function(df, feat, title, low_col) {
  df <- df %>% arrange(padj); df$.feat <- df[[feat]]
  ordered <- reorder(df$.feat, df$fisher_p)
  last_sig <- levels(ordered)[max(which(df$padj < 0.05))]
  vline_x  <- which(levels(ordered) == last_sig) + 1
  ggplot(df, aes(x = reorder(.feat, fisher_p), y = fisher_OR, fill = -log10(fisher_p))) +
    geom_point(shape = 21, aes(size = sig + not_sig, stroke = ifelse(padj < 0.05, 1, 0))) +
    geom_vline(xintercept = vline_x, linetype = "dashed", color = "red") +
    scale_fill_gradient(high = "lightyellow", low = low_col, name = "-log10(p-value)", trans = "reverse") +
    scale_size_continuous(name = "Proteins", range = c(0.5, 3)) +
    labs(title = title, x = NULL, y = "Odds ratio") +
    theme_bw() +
    theme(legend.position = c(0.6, 1), legend.direction = "horizontal",
          legend.background = element_blank(), legend.key = element_rect(fill = "transparent"),
          legend.key.size = unit(2, "mm"), legend.title = element_text(size = 6),
          legend.text = element_text(size = 6), legend.justification = "top", axis.title = element_text(size = 6),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6), text = element_text(size = 6), panel.grid = element_blank()) +
    guides(fill = guide_colorbar(order = 1), size = guide_legend(order = 2),
           stroke = guide_legend(order = 3, override.aes = list(size = 4)))
}

load(here("data", "sig_vs_tissue_cnt_olink.rda"))
pb <- enrichment_dotplot(sig_vs_tissue_cnt_olink, "tissue", "Plasma - Tissue Enrichment", "#C8A2C8")
ggsave(file.path(OUT, "S9b_plasma_tissue_enrichment.pdf"), pb, width = 16.5, height = 6.5, units = "cm", dpi = 600)

load(here("data", "sig_vs_celltype_cnt_olink.rda"))
pc <- enrichment_dotplot(sig_vs_celltype_cnt_olink, "celltype", "Plasma - Cell Type Enrichment", "#C8A2C8")
ggsave(file.path(OUT, "S9c_plasma_celltype_enrichment.pdf"), pc, width = 16.5, height = 6.5, units = "cm", dpi = 600)

# source data now built by R/source_data/ (per-figure script)
cat("Supplementary Figure 9 (plasma tissue/cell-type) written to figures/supplementary/S9_*.pdf\n")
