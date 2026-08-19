# =============================================================================
# Source Data — Supplementary Figure 2 (saliva tissue / cell-type enrichment).
#
# Input:  data/top_tissue_per_protein_saliva.rda, data/sig_vs_tissue_cnt_saliva.rda,
#         data/sig_vs_celltype_cnt_saliva.rda
# Output: source_data/SourceData_SupplementaryFigure2.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)

# a: tissue-binned per-protein table (|beta| + first regulation direction)
load(here("data", "top_tissue_per_protein_saliva.rda"))
a_tissue_binned <- top_tissue_per_protein_saliva %>%
  dplyr::mutate(b_high_numeric = abs(as.numeric(b.high))) %>%
  dplyr::filter(!grepl("^NA.", Assay), !is.na(b_high_dir)) %>% as.data.frame()

# b, c: Fisher tissue / cell-type enrichment count tables (as produced by analysis 05/08/09)
load(here("data", "sig_vs_tissue_cnt_saliva.rda"))
load(here("data", "sig_vs_celltype_cnt_saliva.rda"))

writexl::write_xlsx(list(a_tissue_binned = a_tissue_binned,
                         b_tissue_enrichment = as.data.frame(sig_vs_tissue_cnt_saliva),
                         c_celltype_enrichment = as.data.frame(sig_vs_celltype_cnt_saliva)),
                    here("source_data", "SourceData_SupplementaryFigure2.xlsx"))
cat("Source Data Supplementary Figure 2 written\n")
