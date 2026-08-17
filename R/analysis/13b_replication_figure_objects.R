# =============================================================================
# Analysis 13b — Replication cohort: assemble the Fig. 5 / ED 2-3 figure objects.
#
# Derives the discovery-vs-replication comparison tables and the effect-size
# cluster / enrichment objects that Fig. 5, ED 2 and ED 3 load, from the fitted
# replication models (analyses 11/12/13) and the discovery Olink models (04).
#
# Input:  data/res.olink.linear.rda (discovery),
#         data/res.validation.linear.pilot.rda, data/res.validation.three_way.rda,
#         data/combined.pilot.res.linear.rda (replication models),
#         data/validation.exerome.dat.rda, data/validation_prot.label.rda
# Output: data/matches_all.rda, data/betas_long.rda, data/betas_sig.rda,
#         data/labels_df.rda, data/agree_df.rda, data/go_res_sex_time.rda,
#         data/validation_effect_size_clusters.rda,
#         data/heatmap_timegroup_significant_cluster_enrichment.rda
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(clusterProfiler); library(org.Hs.eg.db)
})
source(here("R/package_loading.R")); source(here("R/functions_loading.R"))

load(here("data", "res.olink.linear.rda"))
load(here("data", "res.validation.linear.pilot.rda"))
load(here("data", "res.validation.three_way.rda"))
load(here("data", "combined.pilot.res.linear.rda"))
load(here("data", "validation.exerome.dat.rda"))
load(here("data", "validation_prot.label.rda"))

# ---- discovery vs replication matched betas ---------------------------------
exerome    <- res.olink.linear
validation <- res.validation.linear.pilot %>% dplyr::rename(Assay = outcome)
matches_all <- validation %>%
  dplyr::inner_join(exerome, by = "Assay", suffix = c(".validation", ".exerome")) %>%
  dplyr::mutate(
    sig_exerome = !is.na(fdr.aov) & fdr.aov < 0.05,
    sig_valid   = !is.na(pval.aov.t.factor.adj) & pval.aov.t.factor.adj < 0.05,
    sig_cat = dplyr::case_when(
      sig_exerome &  sig_valid ~ "Both",
      sig_exerome & !sig_valid ~ "Discovery only",
      !sig_exerome & sig_valid ~ "Replication only",
      TRUE ~ "Neither"),
    sig_cat = factor(sig_cat, levels = c("Both", "Discovery only", "Replication only", "Neither")))
save(matches_all, file = here("data", "matches_all.rda"))

# long form: discovery (x) vs replication (y) beta per time point
betas_long <- dplyr::bind_rows(
  matches_all %>% dplyr::transmute(Assay, Time_Point = "0h",   Beta_Exerome = beta.exposure0.exerome,   Beta_Validation = beta.exposure0.validation,   sig_cat),
  matches_all %>% dplyr::transmute(Assay, Time_Point = "0.5h", Beta_Exerome = beta.exposure0.5.exerome, Beta_Validation = beta.exposure0.5.validation, sig_cat),
  matches_all %>% dplyr::transmute(Assay, Time_Point = "24h",  Beta_Exerome = beta.exposure24.exerome,  Beta_Validation = beta.exposure24.validation,  sig_cat)) %>%
  dplyr::filter(!is.na(Beta_Exerome) & !is.na(Beta_Validation)) %>%
  dplyr::mutate(Time_Point = factor(Time_Point, levels = c("0h", "0.5h", "24h")))
save(betas_long, file = here("data", "betas_long.rda"))

betas_sig <- betas_long %>% dplyr::filter(sig_cat != "Neither") %>% droplevels()
save(betas_sig, file = here("data", "betas_sig.rda"))

# label frame: top 5 up/down concordant per time + opposing extremes + manual set
top_labels <- betas_sig %>% dplyr::filter(sig_cat == "Both") %>%
  dplyr::group_by(Time_Point) %>% dplyr::mutate(effect_mean = (Beta_Exerome + Beta_Validation) / 2) %>%
  dplyr::arrange(desc(effect_mean)) %>% dplyr::slice_head(n = 5) %>%
  dplyr::bind_rows(betas_sig %>% dplyr::filter(sig_cat == "Both") %>%
    dplyr::group_by(Time_Point) %>% dplyr::mutate(effect_mean = (Beta_Exerome + Beta_Validation) / 2) %>%
    dplyr::arrange(effect_mean) %>% dplyr::slice_head(n = 5)) %>%
  dplyr::distinct(Time_Point, Assay, .keep_all = TRUE) %>% dplyr::ungroup()
opp_labels <- betas_sig %>%
  dplyr::filter((Beta_Exerome < 0 & Beta_Validation > 0) | (Beta_Exerome > 0 & Beta_Validation < 0)) %>%
  dplyr::mutate(quadrant = dplyr::case_when(Beta_Exerome < 0 & Beta_Validation > 0 ~ "top_left", TRUE ~ "bottom_right"),
                radius = sqrt(Beta_Exerome^2 + Beta_Validation^2)) %>%
  dplyr::group_by(Time_Point, quadrant) %>% dplyr::slice_max(order_by = radius, n = 2, with_ties = FALSE) %>% dplyr::ungroup()
manual_df <- betas_sig %>% dplyr::filter(Assay %in% c("GH1", "AMBN", "NT5C1A", "AGR2", "CDK11A", "KCTD2", "SASH3", "IL6"))
labels_df <- dplyr::bind_rows(
  top_labels %>% dplyr::select(Assay, Time_Point, Beta_Exerome, Beta_Validation, sig_cat),
  opp_labels %>% dplyr::select(Assay, Time_Point, Beta_Exerome, Beta_Validation, sig_cat),
  manual_df  %>% dplyr::select(Assay, Time_Point, Beta_Exerome, Beta_Validation, sig_cat)) %>%
  dplyr::distinct(Time_Point, Assay, .keep_all = TRUE)
save(labels_df, file = here("data", "labels_df.rda"))

# concordant vs discordant proportion per time point
agree_df <- betas_sig %>%
  dplyr::mutate(Agreement = dplyr::case_when(
    (Beta_Exerome > 0 & Beta_Validation > 0) | (Beta_Exerome < 0 & Beta_Validation < 0) ~ "Concordant",
    TRUE ~ "Discordant")) %>%
  dplyr::count(Time_Point, Agreement, name = "Count") %>%
  dplyr::group_by(Time_Point) %>% dplyr::mutate(Percent = 100 * Count / sum(Count)) %>% dplyr::ungroup()
save(agree_df, file = here("data", "agree_df.rda"))

# ---- sex x time interaction GO enrichment (ED 2) ----------------------------
protein_columns <- intersect(names(validation.exerome.dat), validation_prot.label$Assay)
significant_proteins <- res.validation.three_way %>%
  dplyr::filter(pval.aov.sex.time.adj <= 0.05) %>% dplyr::pull(outcome)
significant_protein_columns <- intersect(protein_columns, significant_proteins)
go_res <- clusterProfiler::enrichGO(
  gene = significant_protein_columns, OrgDb = org.Hs.eg.db,
  universe = unique(validation_prot.label$Assay), ont = "ALL",
  keyType = "SYMBOL", pvalueCutoff = 0.05)
go_res_sex_time <- as.data.frame(go_res)
save(go_res_sex_time, file = here("data", "go_res_sex_time.rda"))

# ---- effect-size k-means clusters + per-cluster GO:BP enrichment (ED 3) -----
plot_data_long <- combined.pilot.res.linear %>%
  dplyr::mutate(group = ifelse(group == "U", "A", group)) %>%
  dplyr::filter(!is.na(pval.aov.group.adj)) %>%
  dplyr::mutate(outcome = factor(outcome, levels = unique(outcome))) %>%
  dplyr::select(-dplyr::starts_with("pval.exposure"), dplyr::starts_with("pval.exposure") & dplyr::ends_with(".adj")) %>%
  dplyr::rename_with(~ gsub("\\.adj$", "", .), dplyr::starts_with("pval.exposure")) %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("beta.exposure") | dplyr::starts_with("se.exposure") | dplyr::starts_with("pval.exposure"),
                      names_to = c(".value", "time_point"), names_pattern = "(beta|se|pval)\\.exposure([0-9.]*)") %>%
  dplyr::mutate(time_point = dplyr::recode(time_point, "0" = "0h", "0.5" = "0.5h", "24" = "24h"))
significant_outcomes <- plot_data_long %>%
  dplyr::filter(pval.aov.t.factor.group.adj < 0.05 | pval.aov.group.adj < 0.05 | pval.aov.t.factor.adj < 0.05) %>%
  dplyr::pull(outcome) %>% unique()
heatmap_matrix <- plot_data_long %>%
  dplyr::filter(outcome %in% significant_outcomes) %>%
  dplyr::select(outcome, group, time_point, beta) %>%
  tidyr::pivot_wider(names_from = c(group, time_point), values_from = beta) %>%
  tibble::column_to_rownames("outcome") %>% as.matrix()

set.seed(123)
km <- kmeans(scale(heatmap_matrix), centers = 10)
cluster_order <- factor(km$cluster, levels = 1:10)
validation_effect_size_clusters <- data.frame(Assay = names(cluster_order),
                                               Cluster = as.numeric(cluster_order))
save(validation_effect_size_clusters, file = here("data", "validation_effect_size_clusters.rda"))

# relabel clusters to the published order, then GO:BP enrichment per cluster
level_mapping <- setNames(1:10, c(8, 7, 2, 9, 10, 4, 1, 6, 5, 3))
cluster_lists <- split(names(cluster_order),
                       factor(level_mapping[as.character(as.integer(cluster_order))], levels = 1:10))
background_list <- unique(validation_prot.label$Assay)
run_go_enrichment <- function(genes, cluster_num) {
  go <- clusterProfiler::enrichGO(gene = genes, OrgDb = org.Hs.eg.db, universe = background_list,
                                  pAdjustMethod = "fdr", minGSSize = 2, maxGSSize = 500,
                                  ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05)
  if (is.null(go) || nrow(as.data.frame(go)) == 0) return(NULL)
  go <- clusterProfiler::simplify(go, cutoff = 0.7, by = "p.adjust", select_fun = min)
  as.data.frame(go) %>% dplyr::mutate(cluster = cluster_num, source = "GO:BP")
}
heatmap_timegroup_significant_cluster_enrichment <- dplyr::bind_rows(
  lapply(seq_along(cluster_lists), function(i) run_go_enrichment(cluster_lists[[i]], i)))
save(heatmap_timegroup_significant_cluster_enrichment,
     file = here("data", "heatmap_timegroup_significant_cluster_enrichment.rda"))

cat("Analysis 13b done: wrote matches_all, betas_long/sig, labels_df, agree_df, go_res_sex_time,\n",
    "  validation_effect_size_clusters, heatmap_timegroup_significant_cluster_enrichment.\n")
