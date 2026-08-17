# =============================================================================
# Supplementary Figure 3 — Temporally resolved exercise-induced urinary proteome dynamics (related to Fig. 2).
#
#   a  per-cluster temporal profiles (5 clusters), 5 lowest-q proteins labelled
#   b  top 10 enriched ORA terms across the regulated urinary proteins
#   c  number of significant ORA terms per cluster
#   d  HPA disease-category enrichment per cluster (metabolic disorders)
#
# Input:  data/res.urine.linear.rda, data/exerome.dat.urine.rda, data/prot.label.urine.rda,
#         data/res.enrich.urine.rda, data/res.enrich.urine.cluster.rda,
#         data/top_tissue_per_protein_urine.rda, data/data_hpa_categorized.rda
# Output: figures/supplementary/S3_*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggpubr); library(reshape2); library(scales)
})
source(here("R", "package_loading.R")); source(here("R", "functions_loading.R"))
OUT <- here("figures", "supplementary"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cluster_color_palette <- c("0"="#B0B0B0","1"="#1F78B4","2"="#33A02C","3"="#E31A1C",
                           "4"="#FDBF6F","5"="#FF7F00")

load(here("data", "res.urine.linear.rda"))
load(here("data", "exerome.dat.urine.rda")); exerome.dat <- exerome.dat.urine
load(here("data", "prot.label.urine.rda"));  prot.label <- prot.label.urine

cl.label.exerome <- tibble::tibble(
  cluster = 1:5,
  label = c("Decrease, <baseline at 24h", "Decrease, >baseline at 24h",
            "Increase, <baseline at 24h", "Increase, >baseline at 24h", "Stochastic pattern"),
  cluster.paper = paste0("C", 1:5))

# ---- a: per-cluster temporal profiles ---------------------------------------
SHEETS <- list()
for (k in 1:5) {
  pk <- cluster_profile_lcms(res.urine.linear, exerome.dat, prot.label, cl.label.exerome, cluster_number = k) +
    coord_cartesian(ylim = c(-2, 2.5))
  ggsave(file.path(OUT, sprintf("S3%s_urine_cluster_%d.pdf", letters[k], k)), pk,
         width = 5.5, height = 5, units = "cm", dpi = 600)
  SHEETS[[sprintf("%s_cluster%d", letters[k], k)]] <- pk$data
}

# ---- b: top 10 ORA terms across all regulated proteins ----------------------
load(here("data", "res.enrich.urine.rda"))
top_terms <- subset(res.enrich.urine, significant == TRUE) %>% arrange(p_value) %>% head(10)
top_terms$term_name <- str_wrap(top_terms$term_name, width = 30)
pb <- ggplot(top_terms, aes(x = p_value, y = reorder(term_name, -p_value))) +
  geom_bar(stat = "identity", fill = "black") +
  scale_x_log10(labels = scales::scientific_format()) +
  labs(x = "p-value", y = NULL, title = "Urine Top Enrichment Terms") +
  theme_bw() +
  theme(text = element_text(size = 6), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT, "S3f_urine_overall_enrichment.pdf"), pb, width = 5.5, height = 5, units = "cm")

# ---- c: number of significant ORA terms per cluster -------------------------
load(here("data", "res.enrich.urine.cluster.rda"))
enrichment_counts <- table(res.enrich.urine.cluster$cluster, res.enrich.urine.cluster$result.significant)
enrichment_df <- as.data.frame.matrix(enrichment_counts)
enrichment_df$Cluster <- rownames(enrichment_df)
colnames(enrichment_df) <- c("Num. of signif_enrich", "Cluster")
enrichment_df_melt <- reshape2::melt(enrichment_df, id.vars = "Cluster")
pc <- ggplot(enrichment_df_melt, aes(x = Cluster, y = value)) +
  geom_bar(stat = "identity", fill = "black") +
  labs(x = "Cluster", y = "Significant Enrichment Terms") +
  theme_bw() +
  theme(axis.text = element_text(size = 6), axis.title = element_text(size = 6),
        plot.title = element_text(size = 6), panel.grid = element_blank())
ggsave(file.path(OUT, "S3g_urine_enrichment_per_cluster.pdf"), pc, width = 5.5, height = 5, units = "cm")

# ---- h: selected GO:BP / Reactome terms uniquely enriched per cluster --------
filtered_df <- res.enrich.urine.cluster %>% filter(!str_detect(result.term_name, "root term"))
specific_terms <- c(
  "Signal Transduction", "RHO GTPase cycle", "Signaling by Insulin receptor",
  "Formation of the cornified envelope", "Keratinization", "Elastic fibre formation",
  "Iron uptake and transport", "Platelet activation, signaling and aggregation",
  "Fcgamma receptor (FCGR) dependent phagocytosis", "Neurotrophin signaling pathway")
wrap_text <- function(x, width = 30) sapply(x, function(y) paste(strwrap(y, width = width), collapse = "\n"))
plot_data <- filtered_df %>% dplyr::filter(result.term_name %in% specific_terms)
plot_data$result.term_name <- wrap_text(plot_data$result.term_name)
ph <- ggplot(plot_data, aes(x = reorder(result.term_name, -result.p_value),
                            y = -log10(result.p_value), fill = as.factor(cluster))) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = NULL, y = "-log10(adj.p-val)") +
  scale_fill_manual(values = cluster_color_palette) +
  theme_pubr() +
  theme(plot.title = element_text(size = 6),
        axis.text.x = element_text(size = 6, angle = 70, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 6), text = element_text(size = 6),
        legend.position = "none", axis.ticks = element_blank(), panel.grid = element_blank()) +
  facet_wrap(~ cluster.paper, scales = "free_x", nrow = 1)
ggsave(file.path(OUT, "S3h_urine_selected_ora_terms.pdf"), ph, width = 11, height = 5, units = "cm", dpi = 600)

# ---- i: incident-disease-risk (HR>1) enrichment per cluster -----------------
# uses the original disease-enrichment plot (plot_disease_enrichment)
source(here("R", "supplementary", "helpers_supp_plots.R"))
inc <- read.csv(here("doc", "supplemental_tables", "disease_enrichment_increasing_urine.csv"))
pi <- plot_disease_enrichment(inc, cluster_color_palette,
                              title = "Disease Enrichment Across Clusters (HR>1)")
ggsave(file.path(OUT, "S3i_urine_incident_disease.pdf"), pi, height = 5, width = 11, units = "cm", dpi = 600)

SHEETS[["f_overall_ORA"]]      <- top_terms
SHEETS[["g_cluster_ORA_counts"]] <- enrichment_df_melt
SHEETS[["h_selected_ORA_terms"]] <- plot_data
SHEETS[["i_incident_disease"]]   <- pi$data
# source data now built by R/source_data/ (per-figure script)
cat("Supplementary Figure 3 (urine dynamics) written to figures/supplementary/S3_*.pdf\n")
