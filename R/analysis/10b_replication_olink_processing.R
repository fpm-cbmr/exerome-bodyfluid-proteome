# =============================================================================
# Analysis 10b — Replication cohort (exercise-mode study) Olink processing.
#
# Reads the two validation Olink Explore plates, builds the wide NPX table and
# protein annotation, and produces the group-wise baseline (t = -1) z-normalised
# frame used by the replication linear models (analyses 11/12) and Fig. 5 / ED 2-3.
#
# Raw inputs (place under data-raw/):
#   data-raw/validation_exercise_mode/validation_olink_Plate1.parquet
#   data-raw/validation_exercise_mode/validation_olink_Plate2.parquet
#   data-raw/validation_exercise_mode/validation_metadata.xlsx
#   data-raw/validation_clinical_dataa.xlsx
#   data-raw/pilot_clinical_dataa.xlsx
#
# Output: data/validation_npx_data.rda, data/validation_prot.label.rda,
#         data/validation.exerome.dat.rda
# =============================================================================
suppressMessages({
  library(here); library(arrow); library(dplyr); library(tidyr); library(stringr)
  library(tibble); library(readxl)
  library(AnnotationDbi); library(org.Hs.eg.db); library(biomaRt)
  library(PhosR); library(SummarizedExperiment)
})
source(here("R/package_loading.R")); source(here("R/functions_loading.R"))

# read one validation plate -> wide protein x sample NPX (rows = proteins)
read_plate <- function(parquet) {
  df <- arrow::read_parquet(here(parquet), as_data_frame = TRUE) %>%
    dplyr::mutate(SampleID = gsub(" ", "", SampleID)) %>%
    dplyr::filter(str_detect(SampleType, "SAMPLE"), str_detect(AssayType, "assay")) %>%
    dplyr::mutate(subject = str_extract(SampleID, "^EX\\d+")) %>%
    dplyr::filter(!is.na(subject))
  wide <- df %>%
    dplyr::group_by(Assay, SampleID) %>%
    dplyr::summarize(mean_NPX = mean(NPX), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = SampleID, values_from = mean_NPX) %>%
    dplyr::rename_with(~ gsub("_0", "_REST24", .), everything()) %>%
    as.data.frame()
  rownames(wide) <- wide$Assay; wide$Assay <- NULL
  colnames(wide) <- gsub("X", "", colnames(wide))
  wide
}

wide_df1 <- read_plate("data-raw/validation_exercise_mode/validation_olink_Plate1.parquet")
wide_df2 <- read_plate("data-raw/validation_exercise_mode/validation_olink_Plate2.parquet")
if (!all(rownames(wide_df1) %in% rownames(wide_df2)))
  stop("Row names do not match between the two validation plates.")
wide_df2 <- wide_df2[match(rownames(wide_df1), rownames(wide_df2)), ]
widedf   <- cbind(wide_df1, wide_df2)

# ---- metadata + SummarizedExperiment ----------------------------------------
clinical_data <- readxl::read_xlsx(here("data-raw", "validation_clinical_dataa.xlsx"))
validation_metadata <- readxl::read_xlsx(here("data-raw", "validation_exercise_mode", "validation_metadata.xlsx")) %>%
  dplyr::left_join(clinical_data, by = "id") %>%
  dplyr::mutate(group = str_extract(sample_id, "(?<=_)[A-Z]"),
                group = ifelse(group == "U", "A", group))   # U -> A so it is the design-matrix reference
widedf <- widedf %>% dplyr::select(all_of(validation_metadata$sample_id))
se_olink_validation <- PhosR::PhosphoExperiment(assay = widedf, colData = validation_metadata)
save(se_olink_validation, file = here("data", "se_olink_validation.rda"))

# ---- wide per-sample NPX table ----------------------------------------------
validation_npx_data <- as.data.frame(t(as.data.frame(SummarizedExperiment::assay(se_olink_validation))))
validation_npx_data$sample <- rownames(validation_npx_data)
validation_npx_data <- validation_npx_data %>%
  tidyr::separate(sample, into = c("subject", "group", "time_point"), sep = "_") %>%
  dplyr::mutate(label = as.factor(subject), time = time_point) %>%
  dplyr::select(-subject, -time_point)
rownames(validation_npx_data) <- NULL
save(validation_npx_data, file = here("data", "validation_npx_data.rda"))

# ---- protein annotation ------------------------------------------------------
cols <- setdiff(names(validation_npx_data), c("time", "label", "group"))
validation_prot.label <- data.frame(protein = cols, Assay = cols, stringsAsFactors = FALSE)
validation_prot.label$ensembl_ids <- AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = validation_prot.label$Assay, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ens2uni <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "uniprotswissprot"),
                          filters = "ensembl_gene_id",
                          values = unique(stats::na.omit(validation_prot.label$ensembl_ids)), mart = mart)
validation_prot.label <- merge(validation_prot.label, ens2uni, by.x = "ensembl_ids",
                               by.y = "ensembl_gene_id", all.x = TRUE) %>%
  dplyr::distinct(Assay, .keep_all = TRUE)
save(validation_prot.label, file = here("data", "validation_prot.label.rda"))

# ---- group-wise baseline (t = -1) z-normalisation ---------------------------
pilot_clinical <- readxl::read_xlsx(here("data-raw", "pilot_clinical_dataa.xlsx")) %>%
  dplyr::rename(label = id)
npx.data <- validation_npx_data %>%
  dplyr::mutate(time = dplyr::recode(time, "PRE" = -1, "POST" = 0, "REST30" = 0.5, "REST24" = 24))
matching_columns <- intersect(validation_prot.label$Assay, colnames(npx.data))

exerome.norm <- npx.data %>% dplyr::filter(time == -1) %>%
  tidyr::pivot_longer(cols = all_of(matching_columns), names_to = "Assay", values_to = "value") %>%
  dplyr::group_by(group, Assay) %>%
  dplyr::summarize(mean.value = mean(value, na.rm = TRUE), sd.value = sd(value, na.rm = TRUE), .groups = "drop")

exerome.dat <- npx.data %>%
  tidyr::pivot_longer(cols = all_of(matching_columns), names_to = "Assay", values_to = "value") %>%
  dplyr::left_join(exerome.norm, by = c("group", "Assay")) %>%
  dplyr::mutate(normalized_value = (value - mean.value) / sd.value) %>%
  tidyr::pivot_wider(id_cols = c(label, time, group), names_from = Assay, values_from = normalized_value) %>%
  dplyr::left_join(pilot_clinical, by = "label") %>%
  dplyr::mutate(t.factor = factor(time), across(c(age, total_fat_mass), scale),
                group = ifelse(group == "U", "A", group))

validation.exerome.dat <- exerome.dat
validation.exerome.dat$sex <- factor(ifelse(validation.exerome.dat$sex == 1, "M", "F"), levels = c("M", "F"))
save(validation.exerome.dat, file = here("data", "validation.exerome.dat.rda"))

cat("Analysis 10b done: wrote validation_npx_data, validation_prot.label, validation.exerome.dat (",
    nrow(validation.exerome.dat), "samples,", length(matching_columns), "proteins )\n")
