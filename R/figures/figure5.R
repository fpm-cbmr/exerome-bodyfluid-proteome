# =============================================================================
# Figure 5 — Replication of the exercise-induced plasma proteome across age, sex, and exercise modalities.
#
# Panels:
#   a  study schematic (BioRender; no source data)
#   b  Venn: discovery-time vs replication-time vs replication-mode proteins
#   c  proportion of concordant exercise-responsive proteins per time point
#   d  beta(discovery) vs beta(replication) scatter, faceted by time
#   e  selected enrichment terms for concordant exerkines ((+)/(-) x time)
#   f  exercise-mode effect-size heatmap (10 k-means clusters, set.seed(123))
#   g  selected ORA terms per replication cluster
#
# Inputs (produced by the analysis pipeline):
#   data/res.olink.linear.rda                                   <- 04_plasma_lmm_clusters.R          (discovery time effect)
#   data/res.validation.linear.pilot.rda                        <- 13_package_replication_period.R
#   data/combined.pilot.res.linear.rda                          <- 13_package_replication_period.R
#   data/matches_all.rda, betas_sig.rda, labels_df.rda,
#     data/agree_df.rda, data/validation_effect_size_clusters.rda,
#     data/heatmap_timegroup_significant_cluster_enrichment.rda <- 13b_replication_figure_objects.R
#   results/ora/concordant_ORA_all.csv   (committed g:Profiler ORA of concordant exerkines)
#
# Output: figures/figure_5/*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(ggplot2); library(ggpubr); library(ggrepel); library(readr)
  library(writexl); library(ggvenn); library(ComplexHeatmap); library(circlize)
})
source(here("R/package_loading.R")); source(here("R/functions_loading.R")); source(here("R/figure_defaults.R"))
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
rename <- dplyr::rename; count <- dplyr::count; summarise <- dplyr::summarise
group_by <- dplyr::group_by; arrange <- dplyr::arrange; slice <- dplyr::slice
dir.create(here("figures", "figure_5"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

group_colors <- c(A = "#004488", H = "#BB5566", S = "#DDAA33")
time_color_palette <- c(`0h` = "#56B4E9", `0.5h` = "#009E73", `24h` = "#CC79A7")

# ---- panel b: Venn — discovery time vs replication time vs replication mode --
load(here("data", "res.olink.linear.rda"))
load(here("data", "res.validation.linear.pilot.rda"))
exerome_linear    <- res.olink.linear
validation_linear <- res.validation.linear.pilot
sig <- 0.05
exerome_proteins    <- exerome_linear %>% filter(fdr.aov < sig) %>% pull(Assay) %>% unique()
validation_proteins <- validation_linear %>% filter(pval.aov.t.factor.adj < sig) %>% pull(outcome) %>% unique()
validation_group    <- validation_linear %>%
  filter(pval.aov.group.adj < sig & pval.aov.t.factor.adj > sig) %>% pull(outcome) %>% unique()
venn_data <- list(Discovery_Time = exerome_proteins,
                  Replication_Time = validation_proteins,
                  Replication_Exercise_Mode = validation_group)
pb <- ggvenn(venn_data, fill_color = c("#601A4A", "#EE442F", "#63ACBE"),
             show_percentage = FALSE, stroke_size = 0.5, set_name_size = 2, text_size = 2) +
  ggtitle("Overlap of plasma proteins with time or mode effect") +
  theme(plot.title = element_text(size = 6, hjust = 0.5))
ggsave(here("figures", "figure_5", "b_venn.pdf"), pb, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["b_venn_counts"]] <- data.frame(
  Set = c("Discovery_Time", "Replication_Time", "Replication_Exercise_Mode"),
  N = c(length(exerome_proteins), length(validation_proteins), length(validation_group)))

# ---- panel c: proportion concordant per time ---------------------------------
load(here("data", "agree_df.rda"))
agree_df <- agree_df %>% mutate(Agreement = factor(Agreement, levels = c("Discordant", "Concordant")),
                                Time_Point = factor(Time_Point, levels = c("0h", "0.5h", "24h")))
pc <- ggplot(agree_df, aes(Time_Point, Percent, fill = Agreement)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.15) +
  scale_fill_manual(values = c(Concordant = "#FF8C00", Discordant = "grey60")) +
  labs(title = "Proportion of concordant exerkines", x = NULL,
       y = "Proportion of 2,514 Proteins (%)", fill = NULL) +
  theme_pubr(base_size = 6) +
  theme(text = element_text(size = 6), legend.position = "top", legend.key.size = unit(2, "mm")) + theme_strokes
ggsave(here("figures", "figure_5", "c_concordance.pdf"), pc, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["c_concordance"]] <- agree_df

# ---- panel d: beta(discovery) vs beta(replication) scatter -------------------
load(here("data", "betas_sig.rda")); load(here("data", "labels_df.rda"))
lev <- c("0h", "0.5h", "24h")
betas_sig <- betas_sig %>% mutate(Time_Point = factor(as.character(Time_Point), levels = lev))
labels_df <- labels_df %>% mutate(Time_Point = factor(as.character(Time_Point), levels = lev))
sig_cols <- c(Both = "#FF8C00", `Discovery only` = "#601A4A", `Replication only` = "#EE442F")
present <- intersect(levels(betas_sig$sig_cat), names(sig_cols))
lim <- ceiling(1.05 * max(abs(c(betas_sig$Beta_Exerome, betas_sig$Beta_Validation)), na.rm = TRUE))
pd <- ggplot(betas_sig, aes(Beta_Exerome, Beta_Validation, color = sig_cat)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.2) +
  geom_point(size = 1.1, alpha = 0.95) +
  ggrepel::geom_text_repel(data = labels_df, aes(label = Assay, color = sig_cat),
                           size = 2, box.padding = 0.2, max.overlaps = Inf,
                           segment.color = "grey60", show.legend = FALSE, xlim = c(-3, 3), ylim = c(-3, lim)) +
  facet_wrap(~ Time_Point, nrow = 1) +
  coord_cartesian(xlim = c(-3, 6), ylim = c(-3, 3), clip = "on") +
  scale_x_continuous(expand = expansion(mult = 0)) + scale_y_continuous(expand = expansion(mult = 0)) +
  scale_color_manual(values = sig_cols[present], breaks = present, drop = FALSE, na.value = "#444444") +
  labs(x = expression(beta ~ "(Discovery)"), y = expression(beta ~ "(Replication)"), color = NULL) +
  theme_pubr(base_size = 6) +
  theme(text = element_text(size = 6), legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_blank(), legend.key.size = unit(2, "mm"),
        strip.text = element_text(size = 6), panel.spacing.x = unit(0.8, "mm")) + theme_strokes
ggsave(here("figures", "figure_5", "d_beta_scatter.pdf"), pd, width = 16.5, height = 5, units = "cm", dpi = 600)
SHEETS[["d_beta_scatter"]] <- betas_sig %>% select(Assay, Time_Point, Beta_Exerome, Beta_Validation, sig_cat)

# ---- panel e: selected enrichment terms for concordant exerkines -------------
ora_tidy <- read_csv(here("results", "ora", "concordant_ORA_all.csv"), show_col_types = FALSE)
focus <- tibble::tribble(
  ~label,                                   ~pattern,
  "Muscle contraction",                     "(^|\\b)muscle\\s+contraction(\\b|$)",
  "CS/DS degradation",                      "cs/ds\\s+degradation|chondroitin\\s+sulfate/.?dermatan\\s+sulfate\\s+degradation",
  "platelet activation",                    "(^|\\b)platelet\\s+activation(\\b|$)",
  "translation",                            "(^|\\b)translation(\\b|$)",
  "Signal Transduction",                    "(^|\\b)signal\\s+transduction(\\b|$)",
  "Metabolism of lipids",                   "metabolism\\s+of\\s+lipids(\\s+and\\s+lipoproteins)?",
  "integrin-mediated signaling pathway",    "integrin[-\\s]mediated\\s+signaling\\s+pathway",
  "vesicle-mediated transport",             "vesicle[-\\s]mediated\\s+transport",
  "translational initiation",               "translational\\s+initiation",
  "TNFs bind their physiological receptors","tnfs\\s+bind\\s+their\\s+physiological\\s+receptors",
  "Extracellular matrix organization",      "extracellular\\s+matrix\\s+organization",
  "cytokine production",                    "(^|\\b)cytokine\\s+production(\\b|$)",
  "neutrophil activation",                  "neutrophil\\s+activation",
  "bone resorption",                        "bone\\s+resorption",
  "glial cell activation",                  "glial\\s+cell\\s+activation",
  "lipid transport",                        "lipid\\s+transport",
  "amyloid-beta clearance",                 "amyloid[-\\s]?beta\\s+clearance",
  "GPCR ligand binding",                    "gpcr\\s+ligand\\s+binding",
  "synaptic signaling",                     "synaptic\\s+signaling",
  "neuron development",                     "neuron\\s+development",
  "gastric motility",                       "gastric\\s+motility")
norm <- function(x) str_to_lower(str_squish(str_replace_all(x, " ", " ")))
ora_clean <- ora_tidy %>% filter(!is.na(term_name)) %>%
  mutate(term_name_norm = norm(term_name),
         Time_Point = factor(Time_Point, levels = lev),
         Direction = recode(Direction, pos = "Positive", neg = "Negative"),
         dir_short = if_else(Direction == "Positive", "(+)", "(-)"))
matches <- tidyr::crossing(ora_clean, focus) %>%
  filter(str_detect(term_name_norm, regex(pattern, ignore_case = TRUE))) %>%
  mutate(label = factor(label, levels = focus$label)) %>%
  arrange(Time_Point, Direction, label, p_value) %>%
  distinct(Time_Point, Direction, label, .keep_all = TRUE) %>%
  mutate(gene_ratio = intersection_size / query_size, color_score = -log10(p_value),
         term_label = factor(str_wrap(as.character(label), 32),
                             levels = unique(str_wrap(focus$label, 32))))
pe <- ggplot(matches, aes(dir_short, term_label)) +
  geom_point(aes(size = gene_ratio, color = color_score), alpha = 0.9, stroke = 0) +
  facet_grid(. ~ Time_Point, scales = "free_x", space = "free_x", switch = "x") +
  labs(x = NULL, y = NULL, size = "Gene Ratio", color = "-log10(q.val)",
       title = "Selected enrichment terms for concordant exerkines") +
  scale_size_continuous(range = c(1, 4), limits = c(0, NA)) +
  theme_pubr(base_size = 6) +
  theme(text = element_text(size = 6), panel.grid = element_blank(), legend.key.size = unit(2, "mm"),
        legend.position = "right", panel.spacing.x = unit(2, "mm"),
        strip.placement = "outside", strip.text.x = element_text(size = 5)) + theme_strokes
ggsave(here("figures", "figure_5", "e_concordant_terms.pdf"), pe, width = 11, height = 10, units = "cm", dpi = 600)
SHEETS[["e_concordant_terms"]] <- matches %>% select(label, Time_Point, dir_short, gene_ratio, color_score, term_name)

# ---- panel f: exercise-mode heatmap (10 k-means clusters) --------------------
res_combined <- combined.pilot.res.linear <- get(load(here("data", "combined.pilot.res.linear.rda")))
plot_data <- res_combined %>% mutate(group = ifelse(group == "U", "A", group)) %>%
  filter(!is.na(pval.aov.group.adj)) %>% mutate(outcome = factor(outcome, levels = unique(outcome)))
plot_long <- plot_data %>%
  select(outcome, group, pval.aov.group.adj, pval.aov.t.factor.adj, pval.aov.t.factor.group.adj,
         starts_with("beta.exposure")) %>%
  pivot_longer(starts_with("beta.exposure"), names_to = "time_point", names_prefix = "beta.exposure",
               values_to = "beta") %>%
  mutate(time_point = recode(time_point, "0" = "0h", "0.5" = "0.5h", "24" = "24h"))
sig_out <- plot_long %>%
  filter(pval.aov.t.factor.group.adj < 0.05 | pval.aov.group.adj < 0.05 | pval.aov.t.factor.adj < 0.05) %>%
  pull(outcome) %>% unique()
heatmap_mat <- plot_long %>% filter(outcome %in% sig_out) %>%
  select(outcome, group, time_point, beta) %>%
  pivot_wider(names_from = c(group, time_point), values_from = beta) %>%
  column_to_rownames("outcome") %>% as.matrix()
group_labels <- sapply(strsplit(colnames(heatmap_mat), "_"), `[`, 1)
time_labels  <- sapply(strsplit(colnames(heatmap_mat), "_"), `[`, 2)
set.seed(123)
km <- kmeans(scale(heatmap_mat), centers = 10)
cluster_order <- factor(km$cluster, levels = 1:10)
col_fun <- colorRamp2(seq(-max(abs(heatmap_mat), na.rm = TRUE), max(abs(heatmap_mat), na.rm = TRUE), length = 3),
                      c("blue", "white", "red"))
column_ha <- HeatmapAnnotation(Group = group_labels, Time = time_labels, show_annotation_name = FALSE,
  simple_anno_size = unit(2, "mm"), col = list(Group = group_colors, Time = time_color_palette))
ht <- Heatmap(heatmap_mat, name = "Effect size", col = col_fun, border = TRUE,
  show_parent_dend_line = FALSE, cluster_rows = TRUE, row_split = cluster_order, row_title = NULL,
  cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE, top_annotation = column_ha,
  heatmap_legend_param = list(title = "Z-score"))
pdf(here("figures", "figure_5", "f_mode_heatmap.pdf"), height = 3.94, width = 4.33)
draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right", merge_legends = TRUE)
dev.off()
SHEETS[["f_mode_heatmap_clusters"]] <- data.frame(Assay = names(cluster_order), Cluster = as.integer(cluster_order))

# ---- panel g: selected ORA terms per replication cluster ---------------------
load(here("data", "heatmap_timegroup_significant_cluster_enrichment.rda"))
enrich <- heatmap_timegroup_significant_cluster_enrichment %>%
  mutate(Description = ifelse(Description == "insulin secretion involved in cellular response to glucose stimulus",
                              "insulin secretion (glucose stimulus)", Description))
target_terms <- c("regeneration", "positive regulation of cell development",
  "extracellular matrix organization", "positive regulation of osteoblast differentiation",
  "response to oxidative stress", "organelle disassembly", "leukocyte activation", "cytokine production",
  "somatodendritic compartment", "coated vesicle", "synapse organization",
  "insulin secretion (glucose stimulus)", "neutrophil mediated cytotoxicity", "hormone secretion")
matching_terms <- enrich %>% filter(Description %in% target_terms) %>%
  mutate(Description = factor(str_wrap(Description, 33), levels = str_wrap(target_terms, 33)))
pg <- ggplot(matching_terms, aes(factor(cluster), Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_size_continuous(range = c(0.5, 2), breaks = c(10, 30, 50)) +
  labs(title = NULL, x = NULL, y = NULL, size = "Count") +
  theme_bw(base_size = 6) +
  theme(legend.position = "right", legend.key.size = unit(2, "mm"), text = element_text(size = 6))
ggsave(here("figures", "figure_5", "g_cluster_ora.pdf"), pg, width = 5.5, height = 5, units = "cm", dpi = 600)
SHEETS[["g_cluster_ora"]] <- matching_terms %>% select(cluster, Description, Count, p.adjust)

# Source Data for this figure is built separately by
# R/source_data/figure5_source_data.R (writes source_data/SourceData_Figure5.xlsx).
cat("Figure 5 panels written:", paste(names(SHEETS), collapse = ", "), "\n")
