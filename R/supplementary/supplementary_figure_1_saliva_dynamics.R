# =============================================================================
# Supplementary Figure 1 — Temporally resolved exercise-induced salivary proteome dynamics (related to Fig. 2).
#
#   a  per-cluster temporal profiles (7 clusters), 5 lowest-q proteins labelled
#   b  top 10 enriched ORA terms across all 601 exercise-regulated saliva proteins
#   c  number of significant ORA terms per cluster
#   d  enrichment of saliva proteins mapping to metabolic disorders (HPA)
#   e  FDA-approved drug targets per cluster (counts above / % below)
#
# Input:  data/res.saliva.linear.rda, data/exerome.dat.saliva.rda,
#         data/prot.label.saliva.rda, data/res.enrich.saliva.rda,
#         data/res.enrich.saliva.cluster.rda, data/top_tissue_per_protein_saliva.rda,
#         data/data_hpa_categorized.rda
# Output: figures/supplementary/S1_*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggpubr); library(reshape2); library(scales)
})
source(here("R", "package_loading.R")); source(here("R", "functions_loading.R"))
OUT <- here("figures", "supplementary"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cluster_color_palette <- c("0"="#B0B0B0","1"="#1F78B4","2"="#33A02C","3"="#E31A1C",
                           "4"="#FDBF6F","5"="#FF7F00","6"="#6A3D9A","7"="#CAB2D6")

load(here("data", "res.saliva.linear.rda"))
load(here("data", "exerome.dat.saliva.rda")); exerome.dat <- exerome.dat.saliva
load(here("data", "prot.label.saliva.rda"));  prot.label <- prot.label.saliva

cl.label.exerome <- tibble::tibble(
  cluster = 1:7,
  label = c("Increase, elevated at 24h", "Decrease, baseline at 24h",
            "Increase at 3h, elevated at 24h", "Decrease, reduced at 24h",
            "Increase, to baseline at 1-3h", "Sequential and sustained increase",
            "Stochastic patterns"),
  cluster.paper = paste0("C", 1:7))

# ---- a: per-cluster temporal profiles ---------------------------------------
SHEETS <- list()   # accumulate source-data sheets
for (k in 1:7) {
  pk <- cluster_profile_lcms(res.saliva.linear, exerome.dat, prot.label, cl.label.exerome, cluster_number = k) +
    coord_cartesian(ylim = c(-2, 2.5))
  ggsave(file.path(OUT, sprintf("S1%s_saliva_cluster_%d.pdf", letters[k], k)), pk,
         width = 5.5, height = 5, units = "cm", dpi = 600)
  SHEETS[[sprintf("%s_cluster%d", letters[k], k)]] <- pk$data
}

# ---- b: top 10 ORA terms across all 601 regulated proteins ------------------
load(here("data", "res.enrich.saliva.rda"))
top_terms <- subset(res.enrich.saliva, significant == TRUE) %>% arrange(p_value) %>% head(10)
top_terms$term_name <- str_wrap(top_terms$term_name, width = 30)
pb <- ggplot(top_terms, aes(x = p_value, y = reorder(term_name, -p_value))) +
  geom_bar(stat = "identity", fill = "black") +
  scale_x_log10(labels = scales::scientific_format()) +
  labs(x = "p-value", y = NULL, title = "Saliva Top Enrichment Terms") +
  theme_bw() +
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT, "S1h_saliva_overall_enrichment.pdf"), pb, width = 5.5, height = 5, units = "cm")
SHEETS[["h_overall_ORA"]] <- top_terms

# ---- c: number of significant ORA terms per cluster -------------------------
load(here("data", "res.enrich.saliva.cluster.rda"))
enrichment_counts <- table(res.enrich.saliva.cluster$cluster, res.enrich.saliva.cluster$result.significant)
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
ggsave(file.path(OUT, "S1i_saliva_enrichment_per_cluster.pdf"), pc, width = 5.5, height = 5, units = "cm")
SHEETS[["i_cluster_ORA_counts"]] <- enrichment_df_melt

# ---- d/e: HPA disease-category + FDA drug-target enrichment per cluster ------
load(here("data", "top_tissue_per_protein_saliva.rda"))
top_tissue_per_protein <- top_tissue_per_protein_saliva
load(here("data", "data_hpa_categorized.rda"))
data_hpa_categorized <- data_hpa_categorized %>% dplyr::rename("Assay" = "Gene") %>%
  tidyr::separate_rows(Disease_Category, sep = ",\\s*")
top_tissue_per_protein <- top_tissue_per_protein %>% mutate(cluster = ifelse(is.na(cluster), 0, cluster))
combined_data <- left_join(top_tissue_per_protein, data_hpa_categorized, by = "Assay")

# per-cluster Fisher's exact test for each disease category
fisher_results <- list()
for (current_cluster in unique(combined_data$cluster)) {
  for (category in unique(combined_data$Disease_Category)) {
    temp <- combined_data %>%
      mutate(in_category = ifelse(Disease_Category == category, 1, 0),
             in_cluster  = ifelse(cluster == current_cluster, 1, 0))
    ct <- table(temp$in_category, temp$in_cluster)
    key <- paste(current_cluster, category, sep = "_")
    if (all(dim(ct) == c(2, 2))) {
      ft <- fisher.test(ct)
      fisher_results[[key]] <- tibble::tibble(
        Cluster = current_cluster, Disease_Category = category,
        yes_category_yes_cluster = ct[2, 2],
        p_value = ft$p.value, odds_ratio = ft$estimate)
    } else {
      fisher_results[[key]] <- tibble::tibble(
        Cluster = current_cluster, Disease_Category = category,
        yes_category_yes_cluster = NA, p_value = NA, odds_ratio = NA)
    }
  }
}
fisher_results_df <- bind_rows(fisher_results) %>%
  group_by(Cluster) %>%
  mutate(total_proteins_in_cluster = sum(yes_category_yes_cluster, na.rm = TRUE),
         percentage = (yes_category_yes_cluster / total_proteins_in_cluster) * 100) %>%
  ungroup()

# d: metabolic disorders — odds ratio vs -log10(p)
metabolic_data <- fisher_results_df %>% filter(Disease_Category == "Metabolic disorders")
pd <- ggplot(metabolic_data, aes(x = odds_ratio, y = -log10(p_value))) +
  geom_point(aes(fill = factor(Cluster), size = yes_category_yes_cluster),
             shape = 21, color = "black", stroke = 1) +
  scale_fill_manual(values = cluster_color_palette, name = NULL) +
  scale_size_continuous(name = "Counts", range = c(1, 3)) +
  labs(title = "Enrichment for metabolic disorders", x = "Odds Ratio",
       y = expression(-log[10](p~value))) +
  theme_pubr() +
  theme(text = element_text(size = 6), legend.position = c(0.99, 0.01),
        legend.justification = c("right", "bottom"), legend.key.size = unit(2, "mm"),
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none", size = guide_legend(override.aes = list(stroke = 0.5)))
ggsave(file.path(OUT, "S1j_saliva_metabolic_disorders.pdf"), pd, height = 5, width = 5.5, units = "cm", dpi = 600)
SHEETS[["j_metabolic_disorders"]] <- metabolic_data

# e: FDA approved drug targets — counts (above) / % (below), mirrored
fda_data <- fisher_results_df %>% filter(Disease_Category == "FDA approved drug targets") %>%
  mutate(percentage_negative = -percentage)
pe <- ggplot(fda_data) +
  geom_bar(aes(x = factor(Cluster), y = yes_category_yes_cluster, fill = factor(Cluster)),
           stat = "identity", position = "identity", width = 0.7, color = "black") +
  geom_bar(aes(x = factor(Cluster), y = percentage_negative, fill = factor(Cluster)),
           stat = "identity", position = "identity", width = 0.7, alpha = 0.6, color = "black") +
  scale_fill_manual(values = cluster_color_palette, name = NULL) +
  labs(title = "FDA Approved Drug Targets", x = "Cluster", y = "Counts (above) / % (below)") +
  theme_pubr() + theme(text = element_text(size = 6), legend.position = "none",
                       plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = abs, name = "Counts (above) / % (below)")
ggsave(file.path(OUT, "S1k_saliva_fda_drug_targets.pdf"), pe, height = 5, width = 5.5, units = "cm", dpi = 600)
SHEETS[["k_FDA_drug_targets"]] <- fda_data

# source data now built by R/source_data/ (per-figure script)
cat("Supplementary Figure 1 (saliva dynamics) written to figures/supplementary/S1_*.pdf\n")
