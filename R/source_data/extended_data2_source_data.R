# =============================================================================
# Source Data — Extended Data Fig. 2 (replication sex x time). Built from the
# committed replication objects (analyses 10b/12/13/13b)
#
# Input:  data/corrected_threeway_period.rds, data/validation.exerome.dat.rda,
#         data/validation_prot.label.rda, data/go_res_sex_time.rda,
#         results/period_sensitivity/enrichment_increasing_sex_exercisemode_PERIOD.xlsx
# Output: source_data/SourceData_ExtendedData2.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(stringr); library(readxl); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

# ---- a: sex x time heatmap row clusters (ward.D2, k=3, seed 123) --------------
sexc <- readRDS(here("data", "corrected_threeway_period.rds"))
load(here("data", "validation.exerome.dat.rda")); exerome.dat <- validation.exerome.dat
load(here("data", "validation_prot.label.rda"))
sex_time <- sexc %>% dplyr::filter(sex_time_fdr < 0.05) %>% dplyr::pull(protein)
keep <- intersect(intersect(names(exerome.dat), validation_prot.label$Assay), sex_time)
mat <- t(as.matrix(exerome.dat[, keep, drop = FALSE]))
ann <- data.frame(time = factor(exerome.dat$time, levels = c(-1, 0, 0.5, 24)),
                  sex = factor(exerome.dat$sex), group = factor(exerome.dat$group))
mat <- mat[, order(ann$time, ann$sex, ann$group)]
mat_s <- t(scale(t(mat)))
set.seed(123)
row_cl <- cutree(hclust(dist(mat_s), method = "ward.D2"), k = 3)
SHEETS[["a_sextime_row_clusters"]] <- data.frame(protein = rownames(mat_s), row_cluster = as.integer(row_cl))

# ---- b: sex x time ORA (curated terms) ---------------------------------------
load(here("data", "go_res_sex_time.rda"))
selected_terms <- c("protein kinase binding", "GTPase binding", "GTPase regulator activity",
  "protein domain specific binding", "microtubule organizing center", "centrosome",
  "organelle envelope", "mitochondrion")
SHEETS[["b_sextime_ORA"]] <- as.data.frame(go_res_sex_time) %>%
  dplyr::filter(Description %in% selected_terms) %>%
  dplyr::mutate(Count = as.numeric(Count), p.adjust = as.numeric(p.adjust)) %>%
  dplyr::select(Description, Count, p.adjust)

# ---- c: incident-disease enrichment (HR>1), committed Fisher table --------------
SHEETS[["c_disease_enrichment"]] <- readxl::read_xlsx(
  here("results", "period_sensitivity", "enrichment_increasing_sex_exercisemode_PERIOD.xlsx")) %>%
  dplyr::mutate(Log10P_signed = dplyr::if_else(Odds_Ratio > 1, -log10(P_value), log10(P_value))) %>%
  dplyr::select(Disease, Cluster, P_value, Odds_Ratio, padj, Log10P_signed) %>% as.data.frame()

writexl::write_xlsx(SHEETS, here("source_data", "SourceData_ExtendedData2.xlsx"))
cat("Source Data Extended Data 2 written:", paste(names(SHEETS), collapse = ", "), "\n")
