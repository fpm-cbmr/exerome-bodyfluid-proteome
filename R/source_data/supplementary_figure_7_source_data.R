# =============================================================================
# Source Data — Supplementary Figure 7 (plasma Olink dynamics).
#
# Input:  data/res.olink.linear.rda, data/olink.exerome.dat.rda,
#         data/res.enrich.olink.rda, data/res.enrich.olink.cluster.rda,
#         data/top_tissue_per_protein_olink.rda,
#         doc/supplemental_tables/disease_enrichment_increasing_plasma.csv
# Output: source_data/SourceData_SupplementaryFigure7.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(readr); library(writexl) })
source(here("R", "supplementary", "helpers_supp_plots.R"))   # plot_secretome_enrichment
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

# a-f: per-protein per-timepoint mean z-score for each Olink cluster (OlinkID keyed)
cluster_profiles <- function(dat, res, id_col = "OlinkID") {
  res_cl <- res %>% dplyr::filter(!is.na(cluster)) %>% dplyr::select(dplyr::all_of(id_col), cluster)
  ids <- intersect(res_cl %>% dplyr::pull(id_col), colnames(dat))
  dplyr::bind_rows(lapply(ids, function(m) {
    dat %>% dplyr::group_by(time_label) %>%
      dplyr::summarise(zscore = mean(.data[[m]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(!!id_col := m)
  })) %>% dplyr::left_join(res_cl, by = id_col)
}
load(here("data", "res.olink.linear.rda")); load(here("data", "olink.exerome.dat.rda"))
SHEETS[["a_f_cluster_profiles"]] <- as.data.frame(cluster_profiles(olink.exerome.dat, res.olink.linear, "OlinkID"))

# g: top-10 overall ORA
load(here("data", "res.enrich.olink.rda"))
SHEETS[["g_overall_ORA"]] <- subset(res.enrich.olink, significant %in% c(TRUE, "TRUE")) %>%
  dplyr::arrange(p_value) %>% head(10) %>% as.data.frame()

# h: significant ORA terms per cluster
load(here("data", "res.enrich.olink.cluster.rda"))
ec <- as.data.frame.matrix(table(res.enrich.olink.cluster$cluster, res.enrich.olink.cluster$result.significant))
ec$Cluster <- rownames(ec); colnames(ec) <- c("n", "Cluster")
SHEETS[["h_cluster_ORA_counts"]] <- reshape2::melt(ec, id.vars = "Cluster")

# i: secreted-protein enrichment per cluster (Fisher; helper computes from data/ objects)
load(here("data", "top_tissue_per_protein_olink.rda"))
ttp_sec <- top_tissue_per_protein_olink %>%
  dplyr::left_join(res.olink.linear %>% dplyr::select(Assay, cluster) %>% dplyr::distinct(), by = "Assay") %>%
  dplyr::mutate(cluster = ifelse(is.na(cluster), 0, cluster)) %>% dplyr::distinct(Assay, .keep_all = TRUE)
pal <- c("0"="#B0B0B0","1"="#1F78B4","2"="#33A02C","3"="#E31A1C","4"="#FDBF6F","5"="#FF7F00","6"="#6A3D9A")
SHEETS[["i_secreted_enrichment"]] <- as.data.frame(plot_secretome_enrichment(ttp_sec, pal)$data)

# j: incident-disease enrichment (HR>1)
SHEETS[["j_incident_disease"]] <- readr::read_csv(
  here("doc", "supplemental_tables", "disease_enrichment_increasing_plasma.csv"), show_col_types = FALSE) %>% as.data.frame()

writexl::write_xlsx(SHEETS, here("source_data", "SourceData_SupplementaryFigure7.xlsx"))
cat("Source Data Supplementary Figure 7 written:", paste(names(SHEETS), collapse = ", "), "\n")
