# =============================================================================
# Extended Data Fig. 3 (related to main Fig. 6) — Human genetic links between exercise-regulated exerkines and disease risk, and exerkines influenced by exercise modality and physiological traits
#
# This script draws panels a-c: the Secreted-Exerkine <-> "Other Disease" network
# (coloc of secreted exercise proteins with the "Other disease" GWAS parent term;
# trait nodes coloured by organ-system `category`). Panels d-f (the exercise-mode
# / other-traits network) are drawn separately by extended_data3b.R.
#
# Built from the committed colocalization results via
# R/figures/helpers_coloc_network.R (same recipe as Fig. 6, set.seed(123)).
#
# Inputs:
#   data/res.olink.linear.rda                                       <- 04_plasma_lmm_clusters.R
#   data-raw/phewas_pqtl/newest_coloc_results_resubmission.xlsx      (cis-pQTL colocalization results)
#   data-raw/phewas_pqtl/final_pqtl_gwas_with_corrected_parents.csv  (curated GWAS trait parents)
#   data-raw/hpa_24.tsv                                              (HPA secretome annotation)
#
# Output: figures/extended_data_3/*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(stringr); library(readr); library(readxl)
  library(igraph); library(writexl); library(ggsci)
})
source(here("R/figure_defaults.R"))
source(here("R/figures/helpers_coloc_network.R"))
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate; rename <- dplyr::rename
FIG <- here("figures", "extended_data_3"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(here("source_data"), showWarnings = FALSE)
OUT <- file.path(FIG, "network_other_diseases_secreted"); SHEETS <- list()

# ---- inputs (same curated coloc pipeline as Fig. 6) --------------------------
load(here("data", "res.olink.linear.rda"))
coloc_res     <- read_xlsx(here("data-raw", "phewas_pqtl", "newest_coloc_results_resubmission.xlsx"))
parent_traits <- read.csv(here("data-raw", "phewas_pqtl", "final_pqtl_gwas_with_corrected_parents.csv"))
hpa_clean <- read_tsv(here("data-raw", "hpa_24.tsv"), show_col_types = FALSE) %>%
  rename(gene = Gene, secretome_function = `Secretome function`) %>%
  mutate(gene = str_trim(as.character(gene)),
         secretome_function = str_squish(as.character(secretome_function)),
         secretome_function = na_if(secretome_function, ""),
         secretome_function = if_else(str_to_lower(secretome_function) %in%
                                        c("na", "n/a", "not available", "-", "none"),
                                      NA_character_, secretome_function))

# "Other" traits + the d3 category palette (ggsci d3 category20 over sorted categories)
other_traits <- parent_traits %>% filter(grepl("^Other", parent_term, ignore.case = TRUE))
other_combined <- cn_prepare_coloc_traits(other_traits, coloc_res, res.olink.linear, hpa_clean)
categories <- sort(unique(other_combined$category))
category_palette <- setNames(ggsci::pal_d3("category20")(length(categories)), categories)

disease_combined_other <- other_combined %>% filter(parent_term %in% "Other disease")

res3 <- cn_build_trait_network_full(
  df = disease_combined_other, olink_df = res.olink.linear, out_prefix = OUT,
  title_text = "Secreted Exerkine <-> Other Disease Network",
  trait_palette = category_palette, color_col = "category", use_secreted = TRUE,
  min_mod_plot = 10L, alpha_adj = 0.05, run_fisher = TRUE, fisher_col = "category",
  top_n_proteins = 10, top_n_traits = 5)

SHEETS[["a_network_edges"]] <- igraph::as_data_frame(res3$g, "edges")
SHEETS[["b_top_proteins"]]  <- res3$prot_trait_comp
SHEETS[["c_top_traits"]]    <- res3$trait_cluster_comp
if (!is.null(res3$fisher_df)) SHEETS[["abc_module_trait_enrichment"]] <- res3$fisher_df

# source data now built by R/source_data/_source_data.R
cat("Extended Data 3 (a-c) panels written:", paste(names(SHEETS), collapse = ", "), "\n")
