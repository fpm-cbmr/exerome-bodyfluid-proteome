# =============================================================================
# Source Data — Figure 3. Built from committed data/ analysis objects, independent
# of the figure-rendering script. One sheet per panel.
#
# Input:  data/plasma_tsne.rda, data/res.olink.linear.rda, data/olink.exerome.dat.rda,
#         data/plasma_variance_decomposition.rda, data/top_tissue_per_protein_olink.rda,
#         data/cluster_vs_{tissue,cell}_cnt_olink.rda, data/res.enrich.olink.cluster.rda,
#         doc/supplemental_tables/disease_enrichment_decreasing_plasma.csv
# Output: source_data/SourceData_Figure3.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(stringr); library(readr); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()
load(here("data", "res.olink.linear.rda")); load(here("data", "olink.exerome.dat.rda"))

# ---- a: plasma t-SNE ---------------------------------------------------------
load(here("data", "plasma_tsne.rda")); SHEETS[["a_plasma_tsne"]] <- as.data.frame(plasma_tsne)

# ---- d: plasma temporal cluster dynamics -------------------------------------
merged_p <- olink.exerome.dat %>%
  tidyr::pivot_longer(cols = -c(time_label, subject, replicate, gender, t.factor),
                      names_to = "OlinkID", values_to = "value") %>%
  dplyr::left_join(res.olink.linear, by = "OlinkID") %>% dplyr::filter(!is.na(cluster))
SHEETS[["d_plasma_dynamics"]] <- merged_p %>% dplyr::group_by(time_label, cluster) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE),
                   sd_value = ifelse(dplyr::n() > 1, sd(value, na.rm = TRUE), 0), n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(se = ifelse(n > 1, sd_value / sqrt(n), 0),
                ci_lower = mean_value - 1.96 * se, ci_upper = mean_value + 1.96 * se) %>% as.data.frame()

# ---- g, h: tissue / cell-type Fisher log2(OR) matrices -----------------------
enrich_or_table <- function(cnt, row_key) {
  cnt <- cnt %>% dplyr::filter(!is.na(.data[[row_key]]), .data[[row_key]] != "") %>%
    dplyr::distinct(.data[[row_key]], cluster, .keep_all = TRUE)
  keep <- cnt %>% dplyr::group_by(.data[[row_key]]) %>%
    dplyr::summarise(any_sig = any(padj < 0.05, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(any_sig) %>% dplyr::pull(1)
  or <- cnt %>% dplyr::filter(.data[[row_key]] %in% keep) %>%
    dplyr::select(dplyr::all_of(row_key), cluster, fisher_OR) %>%
    tidyr::pivot_wider(names_from = cluster, values_from = fisher_OR) %>% as.data.frame()
  num <- setdiff(names(or), row_key)
  or[num] <- lapply(or[num], function(v) { v[v == 0] <- NA; log2(v) })
  names(or)[names(or) %in% num] <- paste0("C", num); or
}
load(here("data", "cluster_vs_tissue_cnt_olink.rda")); SHEETS[["g_plasma_tissue"]]   <- enrich_or_table(cluster_vs_tissue_cnt_olink, "tissue")
load(here("data", "cluster_vs_cell_cnt_olink.rda"));   SHEETS[["h_plasma_celltype"]] <- enrich_or_table(cluster_vs_cell_cnt_olink, "celltype")

# ---- f: plasma secretome proportion ------------------------------------------
load(here("data", "top_tissue_per_protein_olink.rda"))
SHEETS[["f_plasma_secretome"]] <- top_tissue_per_protein_olink %>% dplyr::select(OlinkID, Secretome.function) %>%
  dplyr::left_join(res.olink.linear %>% dplyr::select(OlinkID, cluster), by = "OlinkID") %>%
  dplyr::mutate(cluster = tidyr::replace_na(cluster, 0), Secretome.function = dplyr::na_if(Secretome.function, "")) %>%
  dplyr::filter(!is.na(Secretome.function)) %>%
  dplyr::group_by(Secretome.function, cluster) %>% dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Secretome.function) %>% dplyr::mutate(total = sum(count), proportion = count / total) %>%
  dplyr::ungroup() %>% as.data.frame()

# ---- i: plasma cluster disease enrichment (HR<1) -----------------------------
SHEETS[["i_plasma_disease"]] <- readr::read_csv(
  here("doc", "supplemental_tables", "disease_enrichment_decreasing_plasma.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(Cluster = as.character(Cluster), Odds_Ratio = as.numeric(Odds_Ratio),
                P_value = as.numeric(P_value), padj = as.numeric(padj),
                Log10P_signed = dplyr::if_else(Odds_Ratio > 1, -log10(P_value), log10(P_value))) %>%
  dplyr::filter(Cluster != "0") %>% dplyr::mutate(Odds_Ratio = as.character(Odds_Ratio)) %>% as.data.frame()

# ---- b: baseline vs residual variance scatter --------------------------------
load(here("data", "plasma_variance_decomposition.rda"))
SHEETS[["b_variance_scatter"]] <- plasma_variance_decomposition %>%
  dplyr::select(OlinkID, Assay, Baseline_Variance, Residual_Variance) %>%
  dplyr::left_join(res.olink.linear %>% dplyr::select(OlinkID, cluster), by = "OlinkID") %>% as.data.frame()

# ---- c: NPX over time for NADK and FOLR3 -------------------------------------
ex_ids <- res.olink.linear %>% dplyr::filter(Assay %in% c("NADK", "FOLR3")) %>% dplyr::select(OlinkID, Assay)
SHEETS[["c_individual_proteins"]] <- olink.exerome.dat %>%
  dplyr::select(time_label, replicate, dplyr::all_of(ex_ids$OlinkID)) %>%
  tidyr::pivot_longer(dplyr::all_of(ex_ids$OlinkID), names_to = "OlinkID", values_to = "NPX") %>%
  dplyr::left_join(ex_ids, by = "OlinkID") %>% as.data.frame()

# ---- e: plasma selected cluster ORA (GO:BP) ----------------------------------
load(here("data", "res.enrich.olink.cluster.rda"))
specific_terms_p <- c("Fat digestion and absorption", "leukocyte migration", "response to hormone",
  "neuron projection development", "neurotransmitter secretion", "response to auditory stimulus",
  "neutrophil degranulation", "Degradation of the extracellular matrix", "tissue remodeling",
  "vesicle-mediated transport", "Signaling by Rho GTPases", "endosomal transport",
  "regulation of proteolysis", "apoptotic signaling pathway", "Cytokine Signaling in Immune system",
  "angiogenesis", "regulation of transmembrane transport", "Muscle contraction")
SHEETS[["e_plasma_cluster_ORA"]] <- res.enrich.olink.cluster %>%
  dplyr::filter(!stringr::str_detect(result.term_name, "root term")) %>%
  dplyr::group_by(cluster, result.intersection_size) %>%
  dplyr::slice_min(result.p_value, n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>%
  dplyr::filter(result.term_name %in% specific_terms_p) %>%
  dplyr::mutate(cluster.paper = paste0("C", cluster)) %>% as.data.frame()

writexl::write_xlsx(SHEETS, here("source_data", "SourceData_Figure3.xlsx"))
cat("Source Data Figure 3 written:", paste(names(SHEETS), collapse = ", "), "\n")
