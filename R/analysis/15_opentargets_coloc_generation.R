# =============================================================================
# Analysis 15 — OpenTargets cis-pQTL <-> GWAS colocalization table (release 25.09).
#
# Builds the cis-pQTL <-> GWAS-trait colocalization table for the exercise-
# regulated plasma proteins from the Open Targets Platform (release 25.09):
# loads the colocalisation, credible-set and study tables, joins them, keeps
# significant GWAS<->pQTL colocalizations (H4 >= 0.8, cis only), maps genes to
# symbols, and subsets to the significant exercise proteins. The result is the
# coloc table taken in later (R/analysis/16_pqtl_coloc_networks.R,
# figure6.R / extended_data3.R).
#
# NOT part of run_all: heavy (~15 min) and needs externally downloaded inputs.
#
# DATA — download before running (Open Targets Platform, release 25.09):
#   Colocalisation : https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.09/output/colocalisation_coloc/
#   Credible sets  : https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.09/output/credible_set/
#   Study          : https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.09/output/study/
#   Disease/HPO    : https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.09/output/disease_hpo/
#   (Download every file in each folder.)
#
# EXPECTED LOCAL LAYOUT (relative to the working directory):
#   coloc_parquet/         all colocalisation_coloc parquet files
#   crediblesets_parquet/  all credible_set parquet files
#   study/*.parquet        the study table (part-00000-...snappy.parquet)
#   disease_hpo.parquet, disease.parquet
#   all_sig_proteins.csv   exercise-regulated proteins (column `Assay`)
#
# OUTPUT:
#   acuteexercise_allsigproteins_coloc_OTrelease2509_H4_0.8.xlsx
# =============================================================================
suppressMessages({
  library(arrow); library(purrr); library(dplyr); library(tidyverse)
  library(stringr); library(org.Hs.eg.db); library(writexl)
})

# ---- Load the Open Targets tables -------------------------------------------
OT_coloc_dataset_2509 <- list.files("coloc_parquet", pattern = "\\.parquet$", full.names = TRUE) %>%
  set_names() %>% map_dfr(~ read_parquet(.x), .id = "source")

OT_credset_dataset_2509 <- list.files("crediblesets_parquet", pattern = "\\.parquet$", full.names = TRUE) %>%
  set_names() %>% map_dfr(~ read_parquet(.x), .id = "source")

OT_studyid_2509 <- read_parquet(
  "part-00000-1e42fdda-f476-428f-8a26-2c8b5c9a3b1e-c000.snappy.parquet")

# ---- Join coloc + credible sets + study info --------------------------------
coloc <- OT_coloc_dataset_2509 %>%
  dplyr::select(leftStudyLocusId, rightStudyLocusId, numberColocalisingVariants,
                h0, h1, h2, h3, h4, colocalisationMethod, betaRatioSignAverage)

credset <- OT_credset_dataset_2509 %>%
  dplyr::select(studyLocusId, studyId, variantId, beta, pValueMantissa, pValueExponent,
                effectAlleleFrequencyFromSource, standardError, sampleSize, studyType, isTransQtl)

studyinfo <- OT_studyid_2509 %>%
  dplyr::select(studyId, geneId, projectId, traitFromSource, traitFromSourceMappedIds, pubmedId, biosampleId)

colocalizations_annotated <- coloc %>%
  left_join(credset, by = c("leftStudyLocusId" = "studyLocusId")) %>%
  dplyr::rename(left_studyId = studyId, left_variantId = variantId, left_beta = beta,
                left_pValueMantissa = pValueMantissa, left_pValueExponent = pValueExponent,
                left_effectAlleleFrequencyFromSource = effectAlleleFrequencyFromSource,
                left_standardError = standardError, left_sampleSize = sampleSize,
                left_studyType = studyType, left_isTransQtl = isTransQtl) %>%
  left_join(credset, by = c("rightStudyLocusId" = "studyLocusId")) %>%
  dplyr::rename(right_studyId = studyId, right_variantId = variantId, right_beta = beta,
                right_pValueMantissa = pValueMantissa, right_pValueExponent = pValueExponent,
                right_effectAlleleFrequencyFromSource = effectAlleleFrequencyFromSource,
                right_standardError = standardError, right_sampleSize = sampleSize,
                right_studyType = studyType, right_isTransQtl = isTransQtl) %>%
  left_join(studyinfo, by = c("left_studyId" = "studyId")) %>%
  dplyr::rename(left_geneId = geneId, left_projectId = projectId, left_trait = traitFromSource,
                left_trait_mapped = traitFromSourceMappedIds, left_pubmedId = pubmedId) %>%
  left_join(studyinfo, by = c("right_studyId" = "studyId")) %>%
  dplyr::rename(right_geneId = geneId, right_projectId = projectId, right_trait = traitFromSource,
                right_pubmedId = pubmedId)

# significant cis-pQTL <-> GWAS colocalizations (H4 >= 0.8, cis only)
colocalizations_annotated_pqtl <- colocalizations_annotated %>%
  filter((left_studyType == "gwas" & right_studyType == "pqtl") |
         (left_studyType == "pqtl" & right_studyType == "gwas")) %>%
  filter(right_isTransQtl != "TRUE") %>%
  filter(h4 >= 0.8) %>%
  mutate(logH4H3 = log(h4 / h3))

# ENSEMBL gene id -> symbol
colocalizations_annotated_pqtl$right_geneSymbol <- mapIds(
  org.Hs.eg.db, keys = colocalizations_annotated_pqtl$right_geneId,
  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# ---- Subset to the exercise-regulated proteins and write the table ----------
all_sig_proteins <- read_csv("all_sig_proteins.csv")

all_sign_proteins_coloc_pqtl <-
  colocalizations_annotated_pqtl[colocalizations_annotated_pqtl$right_geneSymbol %in% all_sig_proteins$Assay, ]

write_xlsx(all_sign_proteins_coloc_pqtl,
           "acuteexercise_allsigproteins_coloc_OTrelease2509_H4_0.8.xlsx")

# ---- Trait-category annotation -----
# Maps trait IDs to their parent therapeutic area via the Open Targets disease
if (FALSE) {
  OT_disont_2509 <- read_parquet("disease_hpo.parquet")
  OT_dis         <- read_parquet("disease.parquet")

  key <- OT_dis %>% select(id, name)
  OT_dis_efo   <- OT_dis %>% filter(str_starts(id, "EFO"))
  OT_dis_mondo <- OT_dis %>% filter(str_starts(id, "MONDO"))

  if (!requireNamespace("gwasrapidd", quietly = TRUE)) BiocManager::install("gwasrapidd")
  library(gwasrapidd)

  all_ids <- OT_studyid_2509 %>% pull(traitFromSourceMappedIds) %>%
    flatten_chr() %>% unique() %>% str_replace_all("_", ":")
  efo_ids_list   <- grep("^EFO:",   all_ids, value = TRUE)
  mondo_ids_list <- grep("^MONDO:", all_ids, value = TRUE)

  efo_traits   <- get_traits(efo_uri  = efo_ids_list)
  mondo_traits <- get_traits(trait_id = mondo_ids_list)
  trait_categories <- bind_rows(
    as_tibble(efo_traits@traits)   %>% select(id = efo_id,   category = parent_trait),
    as_tibble(mondo_traits@traits) %>% select(id = trait_id, category = parent_trait))
}
