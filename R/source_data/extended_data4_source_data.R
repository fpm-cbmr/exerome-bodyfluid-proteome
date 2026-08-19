# =============================================================================
# Source Data — Extended Data Fig. 4 (approved drug targets among secreted coloc
# exerkines). Built from data-raw/approved_drug_targets.csv + HPA secretome +
#
# Output: source_data/SourceData_ExtendedData4.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(stringr)
  library(readr); library(snakecase); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

adt <- readr::read_csv(here("data-raw", "approved_drug_targets.csv"), show_col_types = FALSE)
colnames(adt) <- snakecase::to_snake_case(colnames(adt))
adt <- adt %>% dplyr::rename(drug_ids = drug_i_ds) %>% dplyr::mutate(Assay = sub("_.*", "", uniprot_title))
hpa_clean <- readr::read_tsv(here("data-raw", "hpa_24.tsv"), show_col_types = FALSE) %>%
  dplyr::rename(Assay = Gene, secretome_function = `Secretome function`) %>%
  dplyr::mutate(Assay = str_trim(as.character(Assay)),
                secretome_function = dplyr::na_if(str_squish(as.character(secretome_function)), ""),
                secretome_function = dplyr::if_else(str_to_lower(secretome_function) %in%
                  c("na", "n/a", "not available", "-", "none"), NA_character_, secretome_function))
drug_counts <- function(sig_assays) {
  adt %>% dplyr::filter(Assay %in% sig_assays) %>%
    dplyr::mutate(drug_ids = str_split(drug_ids, "[,;]\\s*")) %>% tidyr::unnest(drug_ids) %>%
    dplyr::mutate(drug_ids = str_trim(drug_ids)) %>% dplyr::filter(drug_ids != "") %>%
    dplyr::distinct(Assay, drug_ids) %>% dplyr::inner_join(hpa_clean, by = "Assay") %>%
    dplyr::filter(!is.na(secretome_function)) %>%
    dplyr::group_by(Assay) %>% dplyr::summarise(n_drugs = dplyr::n_distinct(drug_ids), .groups = "drop") %>%
    dplyr::arrange(desc(n_drugs))
}
load(here("data", "res.olink.linear.rda"))
SHEETS[["a_drug_targets_all"]] <- drug_counts(
  res.olink.linear %>% dplyr::filter(fdr.aov < 0.05) %>% dplyr::pull(Assay) %>% unique())
load(here("data", "res.validation.linear.pilot.rda"))
SHEETS[["b_drug_targets_mode"]] <- drug_counts(
  res.validation.linear.pilot %>% dplyr::filter(pval.aov.group.adj < 0.05) %>% dplyr::pull(outcome) %>% unique())

writexl::write_xlsx(SHEETS, here("source_data", "SourceData_ExtendedData4.xlsx"))
cat("Source Data Extended Data 4 written:", paste(names(SHEETS), collapse = ", "), "\n")
