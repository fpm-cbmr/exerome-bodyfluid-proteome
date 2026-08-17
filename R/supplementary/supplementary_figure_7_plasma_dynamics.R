# =============================================================================
# Supplementary Figure 7 — Temporally resolved exercise-induced plasma proteome dynamics (related to Fig. 3).
#
#   a-f  temporally resolved clusters of the Olink plasma proteome (6 clusters),
#        lowest-BH-adj-P proteins coloured, all others grey  [cluster_profile_olink()]
#   g    top 10 ORA terms across all 2,156 exercise-regulated plasma proteins
#   h    number of significant ORA terms per cluster
#   i    enrichment of plasma proteins annotated as secreted, per cluster
#   j    enrichment of plasma proteins linked to increased incident-disease risk
#        (HR>1), per cluster
#
# Uses the original per-cluster plotter cluster_profile_olink() (R/functions_loading.R):
# baseline-normalized Olink NPX, top proteins each in a distinct colour.
#
# Input:  data/res.olink.linear.rda, data/olink.exerome.dat.rda, data/olink.prot.label.rda,
#         data/res.enrich.olink.rda, data/res.enrich.olink.cluster.rda,
#         data/top_tissue_per_protein_olink.rda,
#         doc/supplemental_tables/disease_enrichment_increasing_plasma.csv
# Output: figures/supplementary/S7_*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggpubr); library(ggrepel); library(reshape2); library(scales); library(ggsci)
})
source(here("R", "package_loading.R")); source(here("R", "functions_loading.R"))
OUT <- here("figures", "supplementary"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cluster_color_palette <- c("0"="#B0B0B0","1"="#1F78B4","2"="#33A02C","3"="#E31A1C",
                           "4"="#FDBF6F","5"="#FF7F00","6"="#6A3D9A")

load(here("data", "res.olink.linear.rda"))
load(here("data", "olink.exerome.dat.rda")); exerome.dat <- olink.exerome.dat
load(here("data", "olink.prot.label.rda"));  prot.label <- olink.prot.label

cl.label.olink <- tibble::tibble(
  cluster = 1:6,
  label = c("Decrease, over correction at 1h", "Increase, to baseline at 3h",
            "Increase, over correction", "Increase, slow decrease to 24h",
            "Increase, elevated at 24h", "Stochastic patterns"),
  cluster.paper = paste0("C", 1:6))

# ---- a-f: per-cluster temporal profiles (original generate_plot) -------------
SHEETS <- list()
for (k in 1:6) {
  pk <- cluster_profile_olink(res.olink.linear, exerome.dat, prot.label, cl.label.olink, cluster_number = k)
  ggsave(file.path(OUT, sprintf("S7%s_plasma_cluster_%d.pdf", letters[k], k)), pk,
         width = 5.5, height = 5, units = "cm", dpi = 600)
  SHEETS[[sprintf("%s_cluster%d", letters[k], k)]] <- pk$data
}

# ---- g: top 10 ORA terms across all 2,156 regulated proteins ----------------
load(here("data", "res.enrich.olink.rda"))
top_terms <- subset(res.enrich.olink, significant %in% c(TRUE, "TRUE")) %>% arrange(p_value) %>% head(10)
top_terms$term_name <- str_wrap(top_terms$term_name, width = 30)
pg <- ggplot(top_terms, aes(x = p_value, y = reorder(term_name, -p_value))) +
  geom_bar(stat = "identity", fill = "black") +
  scale_x_log10(labels = scales::scientific_format()) +
  labs(x = "p-value", y = NULL, title = "Plasma Top Enrichment Terms") +
  theme_bw() +
  theme(text = element_text(size = 6), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT, "S7g_plasma_overall_enrichment.pdf"), pg, width = 5.5, height = 5, units = "cm")

# ---- h: number of significant ORA terms per cluster -------------------------
load(here("data", "res.enrich.olink.cluster.rda"))
ec <- table(res.enrich.olink.cluster$cluster, res.enrich.olink.cluster$result.significant)
ec_df <- as.data.frame.matrix(ec); ec_df$Cluster <- rownames(ec_df)
colnames(ec_df) <- c("n", "Cluster")
ph <- ggplot(reshape2::melt(ec_df, id.vars = "Cluster"), aes(Cluster, value)) +
  geom_bar(stat = "identity", fill = "black") +
  labs(x = "Cluster", y = "Significant Enrichment Terms") +
  theme_bw() + theme(axis.text = element_text(size = 6), axis.title = element_text(size = 6), panel.grid = element_blank())
ggsave(file.path(OUT, "S7h_plasma_enrichment_per_cluster.pdf"), ph, width = 5.5, height = 5, units = "cm")

# ---- i: secreted-protein enrichment per cluster (original secretome Fisher) --
source(here("R", "supplementary", "helpers_supp_plots.R"))
load(here("data", "top_tissue_per_protein_olink.rda"))
ttp_sec <- top_tissue_per_protein_olink %>%
  left_join(res.olink.linear %>% select(Assay, cluster) %>% distinct(), by = "Assay") %>%
  mutate(cluster = ifelse(is.na(cluster), 0, cluster)) %>% distinct(Assay, .keep_all = TRUE)
pi <- plot_secretome_enrichment(ttp_sec, cluster_color_palette)
ggsave(file.path(OUT, "S7i_plasma_secreted_enrichment.pdf"), pi, height = 5, width = 5.5, units = "cm", dpi = 600)

# ---- j: incident-disease-risk (HR>1) enrichment (original disease plot) ------
inc <- read.csv(here("doc", "supplemental_tables", "disease_enrichment_increasing_plasma.csv"))
pj <- plot_disease_enrichment(inc, cluster_color_palette,
                              title = "Disease Enrichment Across Clusters (HR>1)")
ggsave(file.path(OUT, "S7j_plasma_incident_disease.pdf"), pj, height = 5, width = 11, units = "cm", dpi = 600)

SHEETS[["g_overall_ORA"]]        <- top_terms
SHEETS[["h_cluster_ORA_counts"]] <- ph$data
SHEETS[["i_secreted_enrichment"]] <- pi$data
SHEETS[["j_incident_disease"]]    <- pj$data
# source data now built by R/source_data/ (per-figure script)
cat("Supplementary Figure 7 (plasma Olink dynamics) written to figures/supplementary/S7_*.pdf\n")
