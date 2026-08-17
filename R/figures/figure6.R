# =============================================================================
# Figure 6 — Human genetics links secreted exerkines to cardiometabolic health and disease.
#
# Panels:
#   a  schematic (BioRender / logos; no source data)
#   b  network: secreted-exerkine cis-pQTL <-> Open GWAS disease
#   c  top traits: composition of connected protein (exercise) clusters
#   d  top proteins: composition of connected trait categories
#   e  trait enrichment per network module (Fisher log2 OR)
#   f  unique GO:BP enriched terms per module
#
# The network is built here from the committed colocalization results via
# R/figures/helpers_coloc_network.R (Louvain communities + Fruchterman-Reingold
# layout, set.seed(123); per-module GO:BP enrichGO is deterministic).
#
# Inputs:
#   data/res.olink.linear.rda                                       <- 04_plasma_lmm_clusters.R  (Olink betas + clusters)
#   data-raw/phewas_pqtl/newest_coloc_results_resubmission.xlsx      (cis-pQTL colocalization results)
#   data-raw/phewas_pqtl/final_pqtl_gwas_with_corrected_parents.csv  (curated GWAS trait parents)
#   data-raw/hpa_24.tsv                                              (HPA secretome annotation)
#   multiomics_II/figures/network_diseases_secreted_testing_Combined_ORA.xlsx  (per-module GO:BP ORA -> f)
#
# Output: figures/figure_6/*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr); library(readr)
  library(readxl); library(tibble); library(ggplot2); library(ggpubr)
  library(igraph); library(scales); library(writexl)
  library(clusterProfiler); library(org.Hs.eg.db); library(fields); library(ggsci)
})
source(here("R/figure_defaults.R"))
source(here("R/figures/helpers_coloc_network.R"))
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
rename <- dplyr::rename; summarise <- dplyr::summarise; group_by <- dplyr::group_by
arrange <- dplyr::arrange; slice_head <- dplyr::slice_head
FIG <- here("figures", "figure_6"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(here("source_data"), showWarnings = FALSE)
OUT <- file.path(FIG, "network_diseases_secreted"); SHEETS <- list()

# ---- inputs & trait palette --------------------------------------------------
load(here("data", "res.olink.linear.rda"))
coloc_res     <- read_xlsx(here("data-raw", "phewas_pqtl", "newest_coloc_results_resubmission.xlsx"))
parent_traits <- read.csv(here("data-raw", "phewas_pqtl", "final_pqtl_gwas_with_corrected_parents.csv"))
all_traits    <- parent_traits
trait_colors <- c(
  "Digestive system disorder" = "#B7704C", "Cardiovascular disease" = "#b33232",
  "Metabolic disorder" = "#fdb462", "Immune system disorder" = "#ffed6f",
  "Neurological disorder" = "#ffffb3", "Lipid or lipoprotein measurement" = "#b3de69",
  "Inflammatory measurement" = "#ccebc5", "Hematological measurement" = "#8dd3c7",
  "Body measurement" = "#66ccff", "Cardiovascular measurement" = "#80b1d3",
  "Other measurement" = "#006699", "Response to drug" = "#fccde5",
  "Biological process" = "#bebada", "Cancer" = "#bc80bd",
  "Other trait" = "#fb8072", "Other disease" = "#ff3399")
cluster_color_palette <- c("0" = "#B0B0B0", "1" = "#1F78B4", "2" = "#33A02C", "3" = "#E31A1C",
                           "4" = "#FDBF6F", "5" = "#FF7F00", "6" = "#6A3D9A", "7" = "#CAB2D6")

hpa_clean <- read_tsv(here("data-raw", "hpa_24.tsv"), show_col_types = FALSE) %>%
  rename(gene = Gene, secretome_function = `Secretome function`) %>%
  mutate(gene = str_trim(as.character(gene)),
         secretome_function = str_squish(as.character(secretome_function)),
         secretome_function = na_if(secretome_function, ""),
         secretome_function = if_else(str_to_lower(secretome_function) %in%
                                        c("na", "n/a", "not available", "-", "none"),
                                      NA_character_, secretome_function))

disease_terms <- c("Cardiovascular disease", "Metabolic disorder", "Immune system disorder",
                   "Cancer", "Digestive system disorder", "Neurological disorder")
disease_combined <- cn_prepare_coloc_traits(all_traits, coloc_res, res.olink.linear, hpa_clean) %>%
  filter(parent_term %in% disease_terms)


res6 <- cn_build_trait_network_full(
  df = disease_combined, olink_df = res.olink.linear, out_prefix = OUT,
  title_text = "Secreted Exerkine cis pQTL <-> Open GWAS Disease Network",
  trait_palette = trait_colors, color_col = "parent_term", use_secreted = TRUE,
  min_mod_plot = 10L, alpha_adj = 0.05, run_fisher = TRUE, fisher_col = "parent_term",
  top_n_proteins = 10, top_n_traits = 5, cluster_color_palette = cluster_color_palette)

SHEETS[["b_network_edges"]] <- igraph::as_data_frame(res6$g, "edges")
SHEETS[["c_top_traits"]]    <- res6$trait_cluster_comp
SHEETS[["d_top_proteins"]]  <- res6$prot_trait_comp
if (!is.null(res6$fisher_df)) SHEETS[["e_module_trait_enrichment"]] <- res6$fisher_df

# ---- panel f: unique GO:BP enriched terms per module (committed ORA workbook) ----
ora_path <- here("multiomics_II", "figures", "network_diseases_secreted_testing_Combined_ORA.xlsx")
if (!file.exists(ora_path)) {
  stop(
    paste0(
      "Missing required ORA workbook for Figure 6 panel f: ", ora_path,
      "\nRegenerate/provide this file before running Figure 6."
    ),
    call. = FALSE
  )
}
ora_combined <- read_xlsx(ora_path)
selected_df <- ora_combined %>% group_by(Module) %>% arrange(p.adjust, desc(Count)) %>%
  slice_head(n = 2) %>% ungroup() %>%
  mutate(log10p = -log10(p.adjust), Module_num = as.numeric(gsub("^M", "", Module)),
         Module = factor(paste0("M", Module_num), levels = paste0("M", sort(unique(Module_num)))),
         Description = str_wrap(ifelse(nchar(Description) > 60, paste0(substr(Description, 1, 57), "..."), Description), 60)) %>%
  arrange(Module_num, desc(log10p), desc(Count), Description) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))
pf <- ggplot(selected_df, aes(Module, Description, size = Count, color = log10p)) +
  geom_point(alpha = 0.85) + scale_color_gradient(low = "lightblue", high = "darkred") +
  scale_size_continuous(range = c(1, 3), breaks = c(1, 4, 8)) + theme_pubr(base_size = 6) +
  labs(title = "Unique GO:BP Enriched Terms per Module", x = "Module", y = NULL,
       size = "Gene Count", color = expression(-log[10]("adj. p"))) +
  theme(text = element_text(size = 6), legend.key.size = unit(2, "mm"), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6), axis.text.y = element_text(size = 5)) + theme_strokes
ggsave(file.path(FIG, "f_unique_go_per_module.pdf"), pf, width = 10, height = 8, units = "cm")
SHEETS[["f_unique_go"]] <- selected_df %>% select(Module, Description, Count, p.adjust, log10p)

# Source Data for this figure is built separately by
# R/source_data/figure6_source_data.R (writes source_data/SourceData_Figure6.xlsx).
cat("Figure 6 panels written:", paste(names(SHEETS), collapse = ", "), "\n")
