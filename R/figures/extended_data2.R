# =============================================================================
# Extended Data Fig. 2 (associated with main Fig. 5) — Proteins with a sex by time interaction in the Replication cohort
#
# Panels:
#   a  sex x time scaled-NPX heatmap (3 row clusters, ward.D2, set.seed(123))
#   b  ORA of the sex x time proteins (curated GO terms, lollipop)
#   c  incident-disease enrichment for Exercise-Mode / Sex-specific proteins (HR>1)
#
# Inputs (produced by the analysis pipeline):
#   data/corrected_threeway_period.rds                     <- 12_replication_sextime.R          (a; sex_time_fdr)
#   data/validation.exerome.dat.rda, data/validation_prot.label.rda
#                                                          <- 10b_replication_olink_processing.R (a)
#   data/go_res_sex_time.rda                               <- 13b_replication_figure_objects.R  (b)
#   results/period_sensitivity/enrichment_increasing_sex_exercisemode_PERIOD.xlsx
#                                (committed Fisher disease-enrichment table)                     (c)
#
# Output: figures/extended_data_2/*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggpubr); library(ggrepel); library(writexl); library(readxl)
  library(ComplexHeatmap); library(circlize)
})
source(here("R/figure_defaults.R"))
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
FIG <- here("figures", "extended_data_2"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

# ---- panel a: sex x time scaled-NPX heatmap ---------------------------------
sexc <- readRDS(here("data", "corrected_threeway_period.rds"))
load(here("data", "validation.exerome.dat.rda")); exerome.dat <- validation.exerome.dat
load(here("data", "validation_prot.label.rda"))
sex_time <- sexc %>% filter(sex_time_fdr < 0.05) %>% pull(protein)
keep <- intersect(intersect(names(exerome.dat), validation_prot.label$Assay), sex_time)
mat <- t(as.matrix(exerome.dat[, keep, drop = FALSE]))
ann <- data.frame(sex = factor(exerome.dat$sex),
                  time = factor(exerome.dat$time, levels = c(-1, 0, 0.5, 24)),
                  group = factor(exerome.dat$group))
ord <- order(ann$time, ann$sex, ann$group); ann <- ann[ord, ]; mat <- mat[, ord]
mat_s <- t(scale(t(mat)))
col_ha <- HeatmapAnnotation(df = ann,
  col = list(sex = c(M = "lightblue", F = "lightpink"),
             time = c(`-1` = "#E69F00", `0` = "#56B4E9", `0.5` = "#009E73", `24` = "#CC79A7"),
             group = c(A = "#1B9E77", H = "#D95F02", S = "#7570B3")),
  simple_anno_size = unit(2, "mm"), annotation_name_gp = gpar(fontsize = 6))
set.seed(123)
row_dend <- hclust(dist(mat_s), method = "ward.D2"); row_cl <- cutree(row_dend, k = 3)
ht <- Heatmap(mat_s, name = "Z-score", top_annotation = col_ha,
  cluster_columns = FALSE, show_column_names = FALSE, show_row_names = FALSE,
  row_split = row_cl, cluster_row_slices = FALSE,
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  column_title = sprintf("Sex x time proteins (corrected + period, n=%d)", length(keep)),
  column_title_gp = gpar(fontsize = 6, fontface = "bold"),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 5), title_gp = gpar(fontsize = 6)))
pdf(file.path(FIG, "a_sextime_heatmap.pdf"), width = 11 / 2.54, height = 8 / 2.54)
draw(ht); dev.off()
SHEETS[["a_sextime_row_clusters"]] <- data.frame(protein = rownames(mat_s), row_cluster = as.integer(row_cl))

# ---- panel b: ORA - Sex*time (curated lollipop) -----------------------------
load(here("data", "go_res_sex_time.rda"))
go <- as.data.frame(go_res_sex_time)
selected_terms <- c("protein kinase binding", "GTPase binding", "GTPase regulator activity",
  "protein domain specific binding", "microtubule organizing center", "centrosome",
  "organelle envelope", "mitochondrion")
filtered_terms <- go %>% filter(Description %in% selected_terms) %>%
  mutate(Count = as.numeric(Count), p.adjust = as.numeric(p.adjust),
         Description = str_wrap(as.character(Description), width = 22))
pb <- ggplot(filtered_terms, aes(y = reorder(Description, -p.adjust), x = -log10(p.adjust))) +
  geom_segment(aes(yend = reorder(Description, -p.adjust), x = 0, xend = -log10(p.adjust)),
               linewidth = 0.5, color = "grey50") +
  geom_point(aes(size = Count), color = "lightpink") +
  scale_size_continuous(range = c(1, 3), breaks = c(20, 28, 36), limits = c(18, 37)) +
  labs(title = "ORA - Sex*time", x = "-log10(Adj.p.val)", y = NULL, size = "Protein count") +
  theme_pubr(base_size = 6) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 6),
        legend.position = c(0.99, 0.28), legend.justification = c("right", "center"),
        legend.key.size = unit(2, "mm")) + theme_strokes
ggsave(file.path(FIG, "b_sextime_ORA.pdf"), pb, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["b_sextime_ORA"]] <- filtered_terms %>% select(Description, Count, p.adjust)

# ---- panel c: incident-disease enrichment (HR>1) ----------------------------
# plotted from the committed Fisher table (avoids the 330 MB ppa_sumstats file)
enr <- read_xlsx(here("results", "period_sensitivity", "enrichment_increasing_sex_exercisemode_PERIOD.xlsx"))
sorted_diseases <- enr %>% filter(Cluster == "Exercise_Mode") %>% arrange(padj) %>% pull(Disease) %>% unique()
enr <- enr %>% mutate(Disease = factor(Disease, levels = sorted_diseases),
                      Log10P_signed = if_else(Odds_Ratio > 1, -log10(P_value), log10(P_value)),
                      Point_Color = if_else(padj < 0.05, as.character(Cluster), "Non-significant"))
cluster_color_palette <- c(Exercise_Mode = "#63ACBE", Sex_Specific = "lightpink", All_Proteins = "grey90")
top_lab <- enr %>% filter(padj < 0.05) %>% group_by(Cluster) %>% arrange(P_value) %>%
  slice_head(n = 3) %>% ungroup() %>%
  bind_rows(enr %>% filter(padj < 0.05 & Odds_Ratio > 1 & Cluster != "All_Proteins")) %>%
  distinct(Cluster, Disease, .keep_all = TRUE)
pc <- ggplot(enr, aes(Disease, Log10P_signed)) +
  geom_point(aes(color = Point_Color, size = Odds_Ratio), alpha = 0.7) +
  scale_color_manual(values = c(cluster_color_palette, `Non-significant` = "grey90")) +
  scale_size_continuous(range = c(0.5, 3), name = "Odds Ratio") +
  geom_text_repel(data = top_lab, aes(Disease, Log10P_signed, label = Disease, color = Cluster),
                  inherit.aes = FALSE, size = 2, force = 40, box.padding = 0.3, max.overlaps = 50) +
  labs(title = "Incident Disease Enrichment for Age/Exercise_Mode Proteins (HR>1)",
       x = "Disease", y = "Signed -log10(P-value)") +
  theme_pubr(base_size = 6) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "right", legend.key.size = unit(2, "mm"), text = element_text(size = 6)) + theme_strokes
ggsave(file.path(FIG, "c_disease_enrichment_increasing.pdf"), pc, width = 11, height = 6, units = "cm", dpi = 600)
SHEETS[["c_disease_enrichment"]] <- enr %>% select(Disease, Cluster, P_value, Odds_Ratio, padj, Log10P_signed)

# source data now built by R/source_data/_source_data.R
cat("Extended Data 2 panels written:", paste(names(SHEETS), collapse = ", "), "\n")
