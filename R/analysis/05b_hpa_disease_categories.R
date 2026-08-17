# =============================================================================
# Analysis 05b — Build HPA disease-category object used by Supplementary Fig. 1.
#
# Creates a gene-level table with disease labels from HPA "Disease involvement"
# and adds an FDA approved drug-target category from DrugBank-derived input.
#
# Input:  data-raw/hpa_24.tsv, data-raw/approved_drug_targets.csv
# Output: data/data_hpa_categorized.rda
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(stringr) })

hpa <- utils::read.delim(here("data-raw", "hpa_24.tsv"), check.names = FALSE)
raw_hpa <- hpa %>%
  transmute(
    Gene = as.character(.data[["Gene"]]),
    Disease = as.character(.data[["Disease involvement"]])
  ) %>%
  filter(!is.na(Gene), Gene != "", !is.na(Disease), Disease != "") %>%
  separate_rows(Disease, sep = ",\\s*") %>%
  mutate(Disease = str_trim(Disease)) %>%
  filter(Disease != "")

adt <- utils::read.csv(here("data-raw", "approved_drug_targets.csv"), check.names = FALSE)
fda_genes <- adt[["Gene Name"]] %>%
  as.character() %>%
  str_trim() %>%
  .[nzchar(.)] %>%
  unique()

fda_hpa <- tibble(
  Gene = fda_genes,
  Disease = "FDA approved drug targets"
)

hpa_long <- bind_rows(raw_hpa, fda_hpa) %>%
  mutate(
    Disease_Category = case_when(
      Disease == "FDA approved drug targets" ~ "FDA approved drug targets",
      str_detect(Disease, regex("diabetes|obesity|lipodystrophy|metabolic", ignore_case = TRUE)) ~ "Metabolic disorders",
      TRUE ~ Disease
    )
  ) %>%
  distinct(Gene, Disease, Disease_Category)

data_hpa_categorized <- hpa_long %>%
  group_by(Gene) %>%
  summarise(
    Disease = paste(sort(unique(Disease)), collapse = ", "),
    Disease_Category = paste(sort(unique(Disease_Category)), collapse = ", "),
    .groups = "drop"
  )

save(data_hpa_categorized, file = here("data", "data_hpa_categorized.rda"))
cat(sprintf("Analysis 05b done: data_hpa_categorized (%d genes).\n", nrow(data_hpa_categorized)))
