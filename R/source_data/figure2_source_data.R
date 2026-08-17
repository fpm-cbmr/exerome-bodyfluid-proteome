# =============================================================================
# Source Data — Figure 2. Built from committed data/ analysis objects, independent
# of the figure-rendering script. One sheet per panel.
#
# Input:  data/{saliva,urine}_tsne.rda, data/res.{saliva,urine}.linear.rda,
#         data/exerome.dat.{saliva,urine}.rda, data/top_tissue_per_protein_saliva.rda,
#         data/res.enrich.saliva.cluster.rda,
#         data/cluster_vs_{tissue,cell}_cnt_{saliva,urine}.rda,
#         doc/supplemental_tables/disease_enrichment_increasing.csv
# Output: source_data/SourceData_Figure2.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(stringr); library(readr); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

# ---- a, h: t-SNE coordinates -------------------------------------------------
load(here("data", "saliva_tsne.rda")); SHEETS[["a_saliva_tsne"]] <- as.data.frame(saliva_tsne)
load(here("data", "urine_tsne.rda"));  SHEETS[["h_urine_tsne"]]  <- as.data.frame(urine_tsne)

# ---- b, i: temporal cluster dynamics (mean z-score + CI per timepoint/cluster)
cluster_dynamics <- function(exerome_dat, res_linear) {
  exerome_dat %>%
    tidyr::pivot_longer(cols = -c(time_label, subject, replicate, gender, t.factor),
                        names_to = "Assay", values_to = "value") %>%
    dplyr::left_join(res_linear, by = "Assay") %>% dplyr::filter(!is.na(cluster)) %>%
    dplyr::group_by(time_label, cluster) %>%
    dplyr::summarise(mean_value = mean(value, na.rm = TRUE),
                     sd_value = ifelse(dplyr::n() > 1, sd(value, na.rm = TRUE), 0),
                     n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(se = ifelse(n > 1, sd_value / sqrt(n), 0),
                  ci_lower = mean_value - 1.96 * se, ci_upper = mean_value + 1.96 * se)
}
load(here("data", "res.saliva.linear.rda")); load(here("data", "exerome.dat.saliva.rda"))
SHEETS[["b_saliva_dynamics"]] <- as.data.frame(cluster_dynamics(exerome.dat.saliva, res.saliva.linear))
load(here("data", "res.urine.linear.rda"));  load(here("data", "exerome.dat.urine.rda"))
SHEETS[["i_urine_dynamics"]] <- as.data.frame(cluster_dynamics(exerome.dat.urine, res.urine.linear))

# ---- c: saliva secreted-protein proportion by cluster ------------------------
load(here("data", "top_tissue_per_protein_saliva.rda"))
SHEETS[["c_saliva_secretome"]] <- top_tissue_per_protein_saliva %>%
  dplyr::mutate(cluster = tidyr::replace_na(cluster, 0),
                Secretome.function = dplyr::na_if(Secretome.function, "")) %>%
  dplyr::filter(!is.na(Secretome.function)) %>%
  dplyr::group_by(Secretome.function, cluster) %>% dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Secretome.function) %>% dplyr::mutate(total = sum(count), proportion = count / total) %>%
  dplyr::ungroup() %>% as.data.frame()

# ---- d: saliva selected cluster ORA (GO:BP) ---------------------------------
load(here("data", "res.enrich.saliva.cluster.rda"))
specific_terms <- list(
  "1" = c("sensory perception of taste", "skeletal system development", "peptide hormone processing", "collagen metabolic process"),
  "6" = c("response to external biotic stimulus", "regulation of vesicle-mediated transport", "programmed cell death"),
  "5" = c("aminoglycan metabolic process", "negative regulation of developmental growth", "Post-translational protein phosphorylation (REAC)"),
  "3" = c("Golgi-to-ER retrograde transport", "Autophagy"),
  "2" = c("Fatty acid degradation", "Glycolysis / Gluconeogenesis", "Lysine degradation"),
  "4" = c("Formation of the cornified envelope", "Keratinization"))
SHEETS[["d_saliva_cluster_ORA"]] <- res.enrich.saliva.cluster %>%
  dplyr::filter(!stringr::str_detect(result.term_name, "root term")) %>%
  dplyr::filter(paste0(cluster) %in% names(specific_terms)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(result.term_name %in% specific_terms[[as.character(unique(cluster))]]) %>%
  dplyr::ungroup() %>% as.data.frame()

# ---- e, f, j: tissue / cell-type Fisher log2(OR) matrices --------------------
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
  names(or)[names(or) %in% num] <- paste0("C", num)
  or
}
load(here("data", "cluster_vs_tissue_cnt_saliva.rda"))
SHEETS[["e_saliva_tissue"]]   <- enrich_or_table(cluster_vs_tissue_cnt_saliva, "tissue")
load(here("data", "cluster_vs_cell_cnt_saliva.rda"))
SHEETS[["f_saliva_celltype"]] <- enrich_or_table(cluster_vs_cell_cnt_saliva, "celltype")
load(here("data", "cluster_vs_cell_cnt_urine.rda"))
SHEETS[["j_urine_celltype"]]  <- enrich_or_table(cluster_vs_cell_cnt_urine, "celltype")

# ---- g: saliva cluster disease enrichment (HR>1) -----------------------------
SHEETS[["g_saliva_disease"]] <- readr::read_csv(
  here("doc", "supplemental_tables", "disease_enrichment_increasing.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(Cluster = as.character(Cluster), Odds_Ratio = as.numeric(Odds_Ratio),
                P_value = as.numeric(P_value), padj = as.numeric(padj),
                Log10P_signed = dplyr::if_else(Odds_Ratio > 1, -log10(P_value), log10(P_value))) %>%
  dplyr::mutate(Odds_Ratio = as.character(Odds_Ratio)) %>% as.data.frame()

writexl::write_xlsx(SHEETS, here("source_data", "SourceData_Figure2.xlsx"))
cat("Source Data Figure 2 written:", paste(names(SHEETS), collapse = ", "), "\n")
