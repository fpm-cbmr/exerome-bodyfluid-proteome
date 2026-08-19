# =============================================================================
# Source Data — Supplementary Figure 9 (plasma tissue / cell-type enrichment).
#
# Input:  data/res.olink.linear.rda, data/top_tissue_per_protein_olink.rda,
#         data/sig_vs_tissue_cnt_olink.rda, data/sig_vs_celltype_cnt_olink.rda
# Output: source_data/SourceData_SupplementaryFigure9.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)

# a: tissue-binned per-protein table (peak |beta| + cluster + fdr joined from the LMM)
load(here("data", "res.olink.linear.rda"))
beta_cols <- grep("^beta\\.exposure", names(res.olink.linear), value = TRUE)
lmm <- res.olink.linear %>%
  dplyr::mutate(b_high_numeric = do.call(pmax, c(lapply(dplyr::across(dplyr::all_of(beta_cols)), abs), na.rm = TRUE))) %>%
  dplyr::select(Assay, cluster, fdr.aov, b_high_numeric) %>% dplyr::distinct(Assay, .keep_all = TRUE)
load(here("data", "top_tissue_per_protein_olink.rda"))
a_tissue_binned <- top_tissue_per_protein_olink %>%
  dplyr::left_join(lmm, by = "Assay") %>%
  dplyr::mutate(cluster = ifelse(is.na(cluster), 0, cluster)) %>%
  dplyr::filter(!grepl("^NA.", Assay), !is.na(b_high_dir), !is.na(fdr.aov)) %>% as.data.frame()

# b, c: Fisher tissue / cell-type enrichment count tables (from analysis 08)
load(here("data", "sig_vs_tissue_cnt_olink.rda"))
load(here("data", "sig_vs_celltype_cnt_olink.rda"))

writexl::write_xlsx(list(a_tissue_binned = a_tissue_binned,
                         b_tissue_enrichment = as.data.frame(sig_vs_tissue_cnt_olink),
                         c_celltype_enrichment = as.data.frame(sig_vs_celltype_cnt_olink)),
                    here("source_data", "SourceData_SupplementaryFigure9.xlsx"))
cat("Source Data Supplementary Figure 9 written\n")
