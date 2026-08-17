# =============================================================================
# Source Data — Supplementary Figure 4 (urine tissue / cell-type enrichment).
# Built from committed data/ analysis objects, independent of the figure.
#
# Input:  data/top_tissue_per_protein_urine.rda, data/sig_vs_tissue_cnt_urine.rda,
#         data/sig_vs_celltype_cnt_urine.rda
# Output: source_data/SourceData_SupplementaryFigure4.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)

# a: tissue-binned per-protein table (first-regulation direction + |beta| from b.high)
load(here("data", "top_tissue_per_protein_urine.rda"))
a_tissue_binned <- top_tissue_per_protein_urine %>%
  dplyr::mutate(b_high_dir = ifelse(as.numeric(b.high) > 0, "up", "down"),
                b_high_numeric = abs(as.numeric(b.high))) %>%
  dplyr::filter(!grepl("^NA.", Assay), !is.na(b_high_dir)) %>% as.data.frame()

# b, c: Fisher tissue / cell-type enrichment count tables (from analysis 05/08/09)
load(here("data", "sig_vs_tissue_cnt_urine.rda"))
load(here("data", "sig_vs_celltype_cnt_urine.rda"))

writexl::write_xlsx(list(a_tissue_binned = a_tissue_binned,
                         b_tissue_enrichment = as.data.frame(sig_vs_tissue_cnt_urine),
                         c_celltype_enrichment = as.data.frame(sig_vs_celltype_cnt_urine)),
                    here("source_data", "SourceData_SupplementaryFigure4.xlsx"))
cat("Source Data Supplementary Figure 4 written\n")
