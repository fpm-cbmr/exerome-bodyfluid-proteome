# =============================================================================
# Source Data — Supplementary Figure 1 (saliva dynamics).
#
# Input:  data/res.saliva.linear.rda, data/exerome.dat.saliva.rda,
#         data/res.enrich.saliva.rda, data/res.enrich.saliva.cluster.rda,
#         data/top_tissue_per_protein_saliva.rda, data/data_hpa_categorized.rda
# Output: source_data/SourceData_SupplementaryFigure1.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(stringr); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

# per-protein per-timepoint mean z-score for every clustered protein (panels a-g)
cluster_profiles <- function(dat, res, id_col = "Assay") {
  res_cl <- res %>% dplyr::filter(!is.na(cluster)) %>% dplyr::select(dplyr::all_of(id_col), cluster)
  ids <- res_cl %>% dplyr::pull(id_col)
  dplyr::bind_rows(lapply(ids, function(m) {
    dat %>% dplyr::group_by(time_label) %>%
      dplyr::summarise(zscore = mean(.data[[m]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(!!id_col := m)
  })) %>% dplyr::left_join(res_cl, by = id_col)
}
load(here("data", "res.saliva.linear.rda")); load(here("data", "exerome.dat.saliva.rda"))
SHEETS[["a_g_cluster_profiles"]] <- as.data.frame(cluster_profiles(exerome.dat.saliva, res.saliva.linear, "Assay"))

# h: top-10 overall ORA terms
load(here("data", "res.enrich.saliva.rda"))
SHEETS[["h_overall_ORA"]] <- subset(res.enrich.saliva, significant == TRUE) %>%
  dplyr::arrange(p_value) %>% head(10) %>% as.data.frame()

# i: number of significant ORA terms per cluster
load(here("data", "res.enrich.saliva.cluster.rda"))
ec <- as.data.frame.matrix(table(res.enrich.saliva.cluster$cluster, res.enrich.saliva.cluster$result.significant))
ec$Cluster <- rownames(ec); colnames(ec) <- c("Num. of signif_enrich", "Cluster")
SHEETS[["i_cluster_ORA_counts"]] <- reshape2::melt(ec, id.vars = "Cluster")

# j, k: per-cluster HPA disease-category Fisher (metabolic disorders + FDA drug targets)
load(here("data", "top_tissue_per_protein_saliva.rda"))
load(here("data", "data_hpa_categorized.rda"))
hpa <- data_hpa_categorized %>% dplyr::rename(Assay = Gene) %>% tidyr::separate_rows(Disease_Category, sep = ",\\s*")
combined <- top_tissue_per_protein_saliva %>% dplyr::mutate(cluster = ifelse(is.na(cluster), 0, cluster)) %>%
  dplyr::left_join(hpa, by = "Assay")
fisher_results <- list()
for (cl in unique(combined$cluster)) for (cat in unique(combined$Disease_Category)) {
  temp <- combined %>% dplyr::mutate(in_category = ifelse(Disease_Category == cat, 1, 0),
                                     in_cluster = ifelse(cluster == cl, 1, 0))
  ct <- table(temp$in_category, temp$in_cluster); key <- paste(cl, cat, sep = "_")
  if (all(dim(ct) == c(2, 2))) { ft <- fisher.test(ct)
    fisher_results[[key]] <- tibble::tibble(Cluster = cl, Disease_Category = cat,
      yes_category_yes_cluster = ct[2, 2], p_value = ft$p.value, odds_ratio = ft$estimate)
  } else fisher_results[[key]] <- tibble::tibble(Cluster = cl, Disease_Category = cat,
      yes_category_yes_cluster = NA, p_value = NA, odds_ratio = NA)
}
fisher_df <- dplyr::bind_rows(fisher_results) %>% dplyr::group_by(Cluster) %>%
  dplyr::mutate(total_proteins_in_cluster = sum(yes_category_yes_cluster, na.rm = TRUE),
                percentage = (yes_category_yes_cluster / total_proteins_in_cluster) * 100) %>% dplyr::ungroup()
SHEETS[["j_metabolic_disorders"]] <- fisher_df %>% dplyr::filter(Disease_Category == "Metabolic disorders") %>% as.data.frame()
SHEETS[["k_FDA_drug_targets"]] <- fisher_df %>% dplyr::filter(Disease_Category == "FDA approved drug targets") %>%
  dplyr::mutate(percentage_negative = -percentage) %>% as.data.frame()

writexl::write_xlsx(SHEETS, here("source_data", "SourceData_SupplementaryFigure1.xlsx"))
cat("Source Data Supplementary Figure 1 written:", paste(names(SHEETS), collapse = ", "), "\n")
