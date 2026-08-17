# =============================================================================
# Replot main figures from publication SourceData workbooks only.
#
# Inputs (required):
#   source_data/SourceData_Figure{1..6}.xlsx
#
# Output:
#   figures_from_source/figure_{1..6}/*.pdf
#
# This script is for scripts-only public releases where intermediate .rda objects
# are not distributed.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(patchwork)
})

OUT <- here("figures_from_source")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

need <- function(path) {
  if (!file.exists(path)) stop(paste("Missing source-data workbook:", path), call. = FALSE)
  path
}

sheet_read <- function(xlsx, sheet) {
  readxl::read_xlsx(xlsx, sheet = sheet)
}

plot_tsne <- function(df, title, out_pdf) {
  p <- ggplot(df, aes(tSNE1, tSNE2, color = as.factor(participant))) +
    geom_point(size = 1) +
    labs(title = title, x = "t-SNE1", y = "t-SNE2", color = "Participant") +
    theme_pubr(base_size = 6)
  ggsave(out_pdf, p, width = 5.5, height = 5, units = "cm", dpi = 600)
}

plot_dynamics <- function(df, ylab, out_pdf) {
  p <- ggplot(df, aes(time_label, mean_value, color = factor(cluster), group = factor(cluster))) +
    geom_line() +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = factor(cluster)), alpha = 0.2, color = NA) +
    geom_hline(yintercept = 0, color = "black") +
    labs(x = "Time [h]", y = ylab, color = "Cluster", fill = "Cluster") +
    theme_pubr(base_size = 6)
  ggsave(out_pdf, p, width = 5.5, height = 5, units = "cm", dpi = 600)
}

plot_marker_timeseries <- function(df, ylab, out_pdf) {
  p <- ggplot(df, aes(time_label, zscore, color = Assay, group = Assay)) +
    geom_line(alpha = 0.8) +
    geom_point(size = 0.7) +
    labs(x = "Time [h]", y = ylab, color = "Assay") +
    theme_pubr(base_size = 6)
  ggsave(out_pdf, p, width = 5.5, height = 5, units = "cm", dpi = 600)
}

plot_secretome <- function(df, out_pdf) {
  p <- ggplot(df, aes(Secretome.function, proportion, fill = factor(cluster))) +
    geom_col() +
    labs(x = "Secretome Function", y = "Proportion", fill = "Cluster") +
    theme_pubr(base_size = 6) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(out_pdf, p, width = 5.5, height = 5, units = "cm", dpi = 600)
}

plot_disease <- function(df, title, out_pdf) {
  p <- ggplot(df, aes(Disease, Log10P_signed)) +
    geom_point(aes(color = as.factor(Cluster), size = Odds_Ratio), alpha = 0.7) +
    labs(title = title, x = "Disease", y = "Signed -log10(P-value)", color = "Cluster", size = "Odds Ratio") +
    theme_pubr(base_size = 6) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave(out_pdf, p, width = 11, height = 6, units = "cm", dpi = 600)
}

plot_dot_enrichment <- function(df, x_col, y_col, size_col, color_col, title, out_pdf, width = 11, height = 8) {
  p <- ggplot(df, aes(.data[[x_col]], .data[[y_col]])) +
    geom_point(aes(size = .data[[size_col]], color = .data[[color_col]]), alpha = 0.85) +
    labs(title = title, x = x_col, y = y_col, size = size_col, color = color_col) +
    theme_pubr(base_size = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out_pdf, p, width = width, height = height, units = "cm", dpi = 600)
}

plot_cluster_ora <- function(df, out_pdf) {
  if (!"cluster.paper" %in% names(df) && "cluster" %in% names(df)) {
    df <- df %>% mutate(cluster.paper = paste0("C", cluster))
  }
  p <- ggplot(df, aes(reorder(result.term_name, -result.p_value), -log10(result.p_value), fill = as.factor(cluster))) +
    geom_col(width = 0.5) +
    facet_wrap(~ cluster.paper, scales = "free_x", nrow = 1) +
    labs(x = NULL, y = "-log10(adj.p-val)") +
    theme_pubr(base_size = 6) +
    theme(axis.text.x = element_text(size = 6, angle = 70, hjust = 1, vjust = 1),
          legend.position = "none", axis.ticks = element_blank())
  ggsave(out_pdf, p, width = 11, height = 5, units = "cm", dpi = 600)
}

plot_venn_counts <- function(df, out_pdf) {
  p <- ggplot(df, aes(Set, N)) +
    geom_col(width = 0.7, fill = "grey40") +
    labs(x = NULL, y = "Count") +
    theme_pubr(base_size = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out_pdf, p, width = 5.5, height = 5, units = "cm", dpi = 600)
}

plot_upset_like <- function(df, out_pdf) {
  if (!all(c("fluids", "Direction") %in% names(df))) return(invisible(NULL))
  bar_df <- df %>% count(fluids, Direction)
  p <- ggplot(bar_df, aes(fluids, n, fill = Direction)) +
    geom_col(width = 0.7) +
    labs(x = NULL, y = "Proteins", fill = NULL) +
    theme_pubr(base_size = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out_pdf, p, width = 5.5, height = 5, units = "cm", dpi = 600)
}

# Figure 1 --------------------------------------------------------------------
fig1_xlsx <- need(here("source_data", "SourceData_Figure1.xlsx"))
fig1_out <- here(OUT, "figure_1")
dir.create(fig1_out, recursive = TRUE, showWarnings = FALSE)

plot_marker_timeseries(sheet_read(fig1_xlsx, "b_metabolic_markers"), "Metabolic marker z-score", file.path(fig1_out, "b_metabolic_markers.pdf"))
plot_marker_timeseries(sheet_read(fig1_xlsx, "c_cytokines"), "Cytokine z-score", file.path(fig1_out, "c_mds_cytokines.pdf"))
{
  d <- sheet_read(fig1_xlsx, "d_protein_counts")
  p <- ggplot(d, aes(factor(time_point), count, fill = fluid)) +
    geom_col(position = "dodge") +
    facet_wrap(~ direction, nrow = 1) +
    labs(x = "Time point", y = "Count", fill = "Fluid") +
    theme_pubr(base_size = 6)
  ggsave(file.path(fig1_out, "d_significant_proteins.pdf"), p, width = 11, height = 5, units = "cm", dpi = 600)
}
plot_upset_like(sheet_read(fig1_xlsx, "e_upset_shared"), file.path(fig1_out, "e_upset_shared_directionality.pdf"))

# Figure 2 --------------------------------------------------------------------
fig2_xlsx <- need(here("source_data", "SourceData_Figure2.xlsx"))
fig2_out <- here(OUT, "figure_2")
dir.create(fig2_out, recursive = TRUE, showWarnings = FALSE)

plot_tsne(sheet_read(fig2_xlsx, "a_saliva_tsne"), "Saliva Proteome Individuality", file.path(fig2_out, "a_saliva_tsne.pdf"))
plot_tsne(sheet_read(fig2_xlsx, "h_urine_tsne"), "Urine Proteome Individuality", file.path(fig2_out, "h_urine_tsne.pdf"))
plot_dynamics(sheet_read(fig2_xlsx, "b_saliva_dynamics"), "Saliva Mean Z-score", file.path(fig2_out, "b_saliva_dynamics.pdf"))
plot_secretome(sheet_read(fig2_xlsx, "c_saliva_secretome"), file.path(fig2_out, "c_saliva_secretome.pdf"))
plot_cluster_ora(sheet_read(fig2_xlsx, "d_saliva_cluster_ORA"), file.path(fig2_out, "d_saliva_cluster_ORA.pdf"))
plot_disease(sheet_read(fig2_xlsx, "g_saliva_disease"), "Saliva disease enrichment", file.path(fig2_out, "g_saliva_disease.pdf"))
plot_dynamics(sheet_read(fig2_xlsx, "i_urine_dynamics"), "Urine Mean Z-score", file.path(fig2_out, "i_urine_dynamics.pdf"))

# Figure 3 --------------------------------------------------------------------
fig3_xlsx <- need(here("source_data", "SourceData_Figure3.xlsx"))
fig3_out <- here(OUT, "figure_3")
dir.create(fig3_out, recursive = TRUE, showWarnings = FALSE)

plot_tsne(sheet_read(fig3_xlsx, "a_plasma_tsne"), "Plasma Proteome Individuality", file.path(fig3_out, "a_plasma_tsne.pdf"))
{
  d <- sheet_read(fig3_xlsx, "b_variance_scatter")
  p <- ggplot(d, aes(Residual_Variance, Baseline_Variance, color = factor(cluster))) +
    geom_point(size = 0.5, alpha = 0.8) +
    labs(x = "Residual variance", y = "Baseline variance", color = "Cluster") +
    theme_pubr(base_size = 6)
  ggsave(file.path(fig3_out, "b_variance_scatter.pdf"), p, width = 5.5, height = 5, units = "cm", dpi = 600)
}
{
  d <- sheet_read(fig3_xlsx, "c_individual_proteins")
  p <- ggplot(d, aes(time_label, NPX, color = factor(replicate), group = interaction(replicate, Assay))) +
    geom_line(alpha = 0.6) +
    facet_wrap(~ Assay, scales = "free_y") +
    labs(x = "Time [h]", y = "NPX", color = "Replicate") +
    theme_pubr(base_size = 6) +
    theme(legend.position = "none")
  ggsave(file.path(fig3_out, "c_individual_proteins.pdf"), p, width = 5.5, height = 5, units = "cm", dpi = 600)
}
plot_dynamics(sheet_read(fig3_xlsx, "d_plasma_dynamics"), "Plasma Mean Z-score", file.path(fig3_out, "d_plasma_dynamics.pdf"))
plot_cluster_ora(sheet_read(fig3_xlsx, "e_plasma_cluster_ORA"), file.path(fig3_out, "e_plasma_cluster_ORA.pdf"))
plot_secretome(sheet_read(fig3_xlsx, "f_plasma_secretome"), file.path(fig3_out, "f_plasma_secretome.pdf"))
plot_disease(sheet_read(fig3_xlsx, "i_plasma_disease"), "Plasma disease enrichment", file.path(fig3_out, "i_plasma_disease.pdf"))

# Figure 4 --------------------------------------------------------------------
fig4_xlsx <- need(here("source_data", "SourceData_Figure4.xlsx"))
fig4_out <- here(OUT, "figure_4")
dir.create(fig4_out, recursive = TRUE, showWarnings = FALSE)

{
  d <- sheet_read(fig4_xlsx, "b_edge_counts")
  p <- ggplot(d, aes(Timepoint)) +
    geom_line(aes(y = TotalEdges, group = 1)) +
    geom_point(aes(y = TotalEdges)) +
    geom_line(aes(y = CrossFluidEdges, group = 1), linetype = "dashed") +
    geom_point(aes(y = CrossFluidEdges), shape = 21, fill = "white") +
    labs(x = "Timepoint", y = "Edges") +
    theme_pubr(base_size = 6)
  ggsave(file.path(fig4_out, "b_edge_counts.pdf"), p, width = 6, height = 5, units = "cm", dpi = 600)
}
{
  d <- sheet_read(fig4_xlsx, "c_MANSC1")
  p <- ggplot(d, aes(Saliva, Plasma, color = Time_Point)) +
    geom_point(size = 1) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.4) +
    labs(x = "Saliva MANSC1 (z-score)", y = "Plasma MANSC1 (z-score)") +
    theme_pubr(base_size = 6)
  ggsave(file.path(fig4_out, "c_mansc1_example.pdf"), p, width = 8.2, height = 4.5, units = "cm", dpi = 600)
}
{
  d <- sheet_read(fig4_xlsx, "f_rho_pre_vs_0h")
  p <- ggplot(d, aes(Time_Pre, Time_0h)) +
    geom_point(size = 0.8, alpha = 0.75) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.4) +
    facet_wrap(~ Source, nrow = 1) +
    labs(x = expression(rho ~ "at Pre"), y = expression(rho ~ "at 0h")) +
    theme_pubr(base_size = 6)
  ggsave(file.path(fig4_out, "f_rho_pre_vs_0h.pdf"), p, width = 16.5, height = 5.5, units = "cm", dpi = 600)
}

# Figure 5 --------------------------------------------------------------------
fig5_xlsx <- need(here("source_data", "SourceData_Figure5.xlsx"))
fig5_out <- here(OUT, "figure_5")
dir.create(fig5_out, recursive = TRUE, showWarnings = FALSE)

{
  d <- sheet_read(fig5_xlsx, "b_venn_counts")
  plot_venn_counts(d, file.path(fig5_out, "b_venn.pdf"))
}
{
  d <- sheet_read(fig5_xlsx, "c_concordance")
  p <- ggplot(d, aes(Time_Point, Percent, fill = Agreement)) +
    geom_col(width = 0.75) +
    labs(x = NULL, y = "Percent", fill = NULL) +
    theme_pubr(base_size = 6)
  ggsave(file.path(fig5_out, "c_concordance.pdf"), p, width = 5.5, height = 5, units = "cm", dpi = 600)
}
{
  d <- sheet_read(fig5_xlsx, "d_beta_scatter")
  p <- ggplot(d, aes(Beta_Exerome, Beta_Validation, color = sig_cat)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_point(size = 1, alpha = 0.8) +
    facet_wrap(~ Time_Point, nrow = 1) +
    labs(x = expression(beta ~ "(Discovery)"), y = expression(beta ~ "(Replication)"), color = NULL) +
    theme_pubr(base_size = 6)
  ggsave(file.path(fig5_out, "d_beta_scatter.pdf"), p, width = 16.5, height = 5, units = "cm", dpi = 600)
}
plot_dot_enrichment(sheet_read(fig5_xlsx, "e_concordant_terms"), "dir_short", "label", "gene_ratio", "color_score", "Concordant term enrichment", file.path(fig5_out, "e_concordant_terms.pdf"), width = 11, height = 10)
plot_dot_enrichment(sheet_read(fig5_xlsx, "g_cluster_ora"), "cluster", "Description", "Count", "p.adjust", "Replication cluster ORA", file.path(fig5_out, "g_cluster_ora.pdf"), width = 5.5, height = 5)

# Figure 6 --------------------------------------------------------------------
fig6_xlsx <- need(here("source_data", "SourceData_Figure6.xlsx"))
fig6_out <- here(OUT, "figure_6")
dir.create(fig6_out, recursive = TRUE, showWarnings = FALSE)

plot_dot_enrichment(sheet_read(fig6_xlsx, "c_top_traits"), "cluster", "Trait", "n", "total_n", "Top traits by cluster", file.path(fig6_out, "c_top_traits.pdf"), width = 11, height = 8)
plot_dot_enrichment(sheet_read(fig6_xlsx, "d_top_proteins"), "TraitGroup", "Protein", "n", "total_n", "Top proteins by trait", file.path(fig6_out, "d_top_proteins.pdf"), width = 11, height = 8)
plot_dot_enrichment(sheet_read(fig6_xlsx, "e_module_trait_enrichment"), "Module", "TraitGroup", "a", "OR", "Module-trait enrichment", file.path(fig6_out, "e_module_trait_enrichment.pdf"), width = 11, height = 8)
plot_dot_enrichment(sheet_read(fig6_xlsx, "f_unique_go"), "Module", "Description", "Count", "log10p", "Unique GO terms", file.path(fig6_out, "f_unique_go_per_module.pdf"), width = 10, height = 8)

message("Replot from source data complete. Output root: ", OUT)
