# =============================================================================
# Analysis 16 — Secreted-exerkine cis-pQTL <-> trait/disease COLOCALIZATION tables
#               (inputs to the Fig 6 / ED 3 genetic networks).
#
# Builds the protein<->trait tables for the Fig 6 / ED 3 genetic networks from
# precomputed cis-pQTL x GWAS-trait colocalization results.
# The coloc step itself is analysis 15 (R/analysis/15_opentargets_coloc_generation.R),
#
# Inputs:
#   data-raw/phewas_pqtl/newest_coloc_results_resubmission.xlsx   (coloc results = output of analysis 15)
#   data-raw/phewas_pqtl/final_pqtl_gwas_with_corrected_parents.csv (curated trait parents)
#   data-raw/hpa_24.tsv                                           (HPA secretome)
#   data/res.olink.linear.rda                                     (Olink clusters + betas)
# Output:
#   data/coloc_trait_disease.rda   (disease parent terms; Fig 6 + ED 3 a-c)
#   data/coloc_trait_other.rda     ("other" traits;       ED 3 d-f uses the mode variant)
#
# The networks themselves (Louvain modules, ORA, layout) are built with the
# shared helper R/figures/helpers_coloc_network.R and rendered in
# figure6.R / extended_data3.R.
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(stringr); library(readr); library(readxl) })
source(here("R/figures/helpers_coloc_network.R"))   # cn_prepare_coloc_traits()
rename <- dplyr::rename; filter <- dplyr::filter; mutate <- dplyr::mutate   # helper masks these

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

disease_terms <- c("Cardiovascular disease", "Metabolic disorder", "Immune system disorder",
                   "Cancer", "Digestive system disorder", "Neurological disorder")

coloc_all <- cn_prepare_coloc_traits(parent_traits, coloc_res, res.olink.linear, hpa_clean)
coloc_trait_disease <- coloc_all %>% filter(parent_term %in% disease_terms)
coloc_trait_other   <- coloc_all %>% filter(grepl("^Other", parent_term, ignore.case = TRUE))
save(coloc_trait_disease, file = here("data", "coloc_trait_disease.rda"))
save(coloc_trait_other,   file = here("data", "coloc_trait_other.rda"))

cat(sprintf("Analysis 15 done: coloc tables | disease: %d rows / %d proteins / %d traits; other: %d rows.\n",
            nrow(coloc_trait_disease), dplyr::n_distinct(coloc_trait_disease$right_geneSymbol),
            dplyr::n_distinct(coloc_trait_disease$left_trait), nrow(coloc_trait_other)))
