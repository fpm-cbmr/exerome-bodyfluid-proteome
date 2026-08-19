# =============================================================================
# Source Data — Supplementary Figure 3 (urine dynamics).
#
# Input:  data/res.urine.linear.rda, data/exerome.dat.urine.rda,
#         data/res.enrich.urine.rda, data/res.enrich.urine.cluster.rda,
#         doc/supplemental_tables/disease_enrichment_increasing_urine.csv
# Output: source_data/SourceData_SupplementaryFigure3.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(stringr); library(readr); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

cluster_profiles <- function(dat, res, id_col = "Assay") {
  res_cl <- res %>% dplyr::filter(!is.na(cluster)) %>% dplyr::select(dplyr::all_of(id_col), cluster)
  ids <- res_cl %>% dplyr::pull(id_col)
  dplyr::bind_rows(lapply(ids, function(m) {
    dat %>% dplyr::group_by(time_label) %>%
      dplyr::summarise(zscore = mean(.data[[m]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(!!id_col := m)
  })) %>% dplyr::left_join(res_cl, by = id_col)
}
load(here("data", "res.urine.linear.rda")); load(here("data", "exerome.dat.urine.rda"))
SHEETS[["a_e_cluster_profiles"]] <- as.data.frame(cluster_profiles(exerome.dat.urine, res.urine.linear, "Assay"))

load(here("data", "res.enrich.urine.rda"))
SHEETS[["f_overall_ORA"]] <- subset(res.enrich.urine, significant == TRUE) %>%
  dplyr::arrange(p_value) %>% head(10) %>% as.data.frame()

load(here("data", "res.enrich.urine.cluster.rda"))
ec <- as.data.frame.matrix(table(res.enrich.urine.cluster$cluster, res.enrich.urine.cluster$result.significant))
ec$Cluster <- rownames(ec); colnames(ec) <- c("Num. of signif_enrich", "Cluster")
SHEETS[["g_cluster_ORA_counts"]] <- reshape2::melt(ec, id.vars = "Cluster")

specific_terms <- c("Signal Transduction", "RHO GTPase cycle", "Signaling by Insulin receptor",
  "Formation of the cornified envelope", "Keratinization", "Elastic fibre formation",
  "Iron uptake and transport", "Platelet activation, signaling and aggregation",
  "Fcgamma receptor (FCGR) dependent phagocytosis", "Neurotrophin signaling pathway")
SHEETS[["h_selected_ORA_terms"]] <- res.enrich.urine.cluster %>%
  dplyr::filter(!stringr::str_detect(result.term_name, "root term")) %>%
  dplyr::filter(result.term_name %in% specific_terms) %>% as.data.frame()

SHEETS[["i_incident_disease"]] <- readr::read_csv(
  here("doc", "supplemental_tables", "disease_enrichment_increasing_urine.csv"), show_col_types = FALSE) %>% as.data.frame()

writexl::write_xlsx(SHEETS, here("source_data", "SourceData_SupplementaryFigure3.xlsx"))
cat("Source Data Supplementary Figure 3 written:", paste(names(SHEETS), collapse = ", "), "\n")
