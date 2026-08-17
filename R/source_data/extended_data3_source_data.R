# =============================================================================
# Source Data — Extended Data Fig. 3 (consolidated a-f).
# Extended Data Fig. 3 is drawn by two scripts (extended_data3.R = panels a-c,
# the secreted-exerkine <-> Other Disease network; extended_data3b.R = panels d-f,
# the exercise-mode <-> Other Traits network). This single source-data script
# rebuilds BOTH networks from the committed analysis inputs and writes ONE file,
# replacing the earlier split SourceData_ExtendedData3.xlsx / _def.xlsx.
#
# Output: source_data/SourceData_ExtendedData3.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(stringr); library(readr)
  library(readxl); library(igraph); library(ggsci)
  library(clusterProfiler); library(org.Hs.eg.db); library(writexl) })
source(here("R/figures/helpers_coloc_network.R"))
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

load(here("data", "res.olink.linear.rda"))
coloc_res     <- readxl::read_xlsx(here("data-raw", "phewas_pqtl", "newest_coloc_results_resubmission.xlsx"))
parent_traits <- read.csv(here("data-raw", "phewas_pqtl", "final_pqtl_gwas_with_corrected_parents.csv"))
hpa_clean <- readr::read_tsv(here("data-raw", "hpa_24.tsv"), show_col_types = FALSE) %>%
  dplyr::rename(gene = Gene, secretome_function = `Secretome function`) %>%
  dplyr::mutate(gene = str_trim(as.character(gene)),
                secretome_function = dplyr::na_if(str_squish(as.character(secretome_function)), ""),
                secretome_function = dplyr::if_else(str_to_lower(secretome_function) %in%
                  c("na", "n/a", "not available", "-", "none"), NA_character_, secretome_function))
other_traits <- parent_traits %>% dplyr::filter(grepl("^Other", parent_term, ignore.case = TRUE))

# ---- panels a-c: secreted exerkine <-> Other Disease network -----------------
other_combined <- cn_prepare_coloc_traits(other_traits, coloc_res, res.olink.linear, hpa_clean)
categories <- sort(unique(other_combined$category))
category_palette <- setNames(ggsci::pal_d3("category20")(length(categories)), categories)
res3 <- cn_build_trait_network_full(
  df = other_combined %>% dplyr::filter(parent_term %in% "Other disease"),
  olink_df = res.olink.linear, out_prefix = file.path(tempdir(), "ed3_abc"),
  title_text = "Secreted Exerkine <-> Other Disease Network",
  trait_palette = category_palette, color_col = "category", use_secreted = TRUE,
  min_mod_plot = 10L, alpha_adj = 0.05, run_fisher = TRUE, fisher_col = "category",
  top_n_proteins = 10, top_n_traits = 5)
SHEETS[["a_network_nodes"]] <- igraph::as_data_frame(res3$g, "vertices")
SHEETS[["a_network_edges"]] <- igraph::as_data_frame(res3$g, "edges")
SHEETS[["b_top_proteins"]]  <- res3$prot_trait_comp
SHEETS[["c_top_traits"]]    <- res3$trait_cluster_comp
if (!is.null(res3$fisher_df)) SHEETS[["abc_module_trait_enrichment"]] <- res3$fisher_df

# ---- panels d-f: exercise-mode <-> Other Traits network ----------------------
load(here("data", "validation_prot.label.rda")); prot.label <- validation_prot.label
load(here("data", "res.validation.linear.pilot.rda"))
grp_clusters <- get(load(here("data", "validation_effect_size_clusters.rda")))
combined.pilot.res.linear <- get(load(here("data", "combined.pilot.res.linear.rda")))
res_grouxtime <- res.validation.linear.pilot %>% dplyr::rename(Assay = outcome) %>%
  dplyr::filter(pval.aov.group.adj < 0.05) %>%
  dplyr::inner_join(grp_clusters, by = "Assay") %>% dplyr::rename(cluster = Cluster)
prepare_coloc_traits_validation <- function(trait_df, coloc_res, res_valid, hpa_clean) {
  trait_map <- trait_df %>% dplyr::rename(reported_traits = dplyr::any_of(c("reported_traits", "reported_trait"))) %>%
    dplyr::select(reported_traits, parent_term, category) %>%
    dplyr::mutate(reported_traits = str_trim(as.character(reported_traits)))
  coloc_df <- coloc_res %>% dplyr::mutate(left_trait = str_trim(as.character(left_trait)),
    right_geneSymbol = str_to_upper(str_trim(as.character(right_geneSymbol)))) %>%
    dplyr::inner_join(trait_map, by = c("left_trait" = "reported_traits"))
  olink_cluster <- res_valid %>% dplyr::mutate(Assay = str_to_upper(str_trim(as.character(Assay)))) %>%
    dplyr::select(Assay, cluster, beta.exposure0) %>% dplyr::distinct(Assay, .keep_all = TRUE)
  coloc_df %>% dplyr::inner_join(olink_cluster, by = c("right_geneSymbol" = "Assay")) %>%
    dplyr::left_join(hpa_clean %>% dplyr::select(gene, secretome_function) %>%
      dplyr::mutate(gene = str_to_upper(str_trim(as.character(gene)))) %>% dplyr::distinct(),
      by = c("right_geneSymbol" = "gene"))
}
other_combined_v <- prepare_coloc_traits_validation(other_traits, coloc_res, res_grouxtime, hpa_clean)
cats_v <- sort(unique(other_combined_v$category))
res3b <- cn_build_trait_network_full(
  df = other_combined_v, olink_df = prot.label, out_prefix = file.path(tempdir(), "ed3_def"),
  title_text = "Secreted Exercise Mode Exerkines <-> Other Traits Network",
  trait_palette = setNames(ggsci::pal_d3("category20")(length(cats_v)), cats_v),
  color_col = "category", use_secreted = TRUE, min_mod_plot = 10L, alpha_adj = 0.05,
  run_fisher = TRUE, fisher_col = "category", top_n_proteins = 10, top_n_traits = 5, emit_cluster_bar = FALSE)
SHEETS[["d_network_nodes"]] <- igraph::as_data_frame(res3b$g, "vertices")
SHEETS[["d_network_edges"]] <- igraph::as_data_frame(res3b$g, "edges")
SHEETS[["e_top_proteins"]]  <- res3b$prot_trait_comp

# panel f: top traits by connected-protein exercise mode
mode_map <- combined.pilot.res.linear %>%
  dplyr::mutate(sig = (pval.exposure0.adj < 0.05) | (pval.exposure0.5.adj < 0.05) | (pval.exposure24.adj < 0.05)) %>%
  dplyr::filter(sig) %>% dplyr::mutate(Assay = str_to_upper(str_trim(as.character(outcome))),
    group = dplyr::recode(group, U = "A")) %>%
  dplyr::group_by(Assay) %>% dplyr::summarise(nmode = dplyr::n_distinct(group), mode1 = dplyr::first(group), .groups = "drop") %>%
  dplyr::mutate(Exercise_Mode = ifelse(nmode > 1, "Multiple", mode1))
SHEETS[["f_top_traits_mode"]] <- res3b$edges_group %>%
  dplyr::filter(Trait %in% res3b$top_traits$name) %>%
  dplyr::left_join(mode_map %>% dplyr::select(Assay, Exercise_Mode), by = c("Protein" = "Assay")) %>%
  dplyr::mutate(Exercise_Mode = ifelse(is.na(Exercise_Mode), "Multiple", Exercise_Mode)) %>%
  dplyr::group_by(Trait, Exercise_Mode) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Trait) %>% dplyr::mutate(total_n = sum(n)) %>% dplyr::ungroup() %>% as.data.frame()

writexl::write_xlsx(SHEETS, here("source_data", "SourceData_ExtendedData3.xlsx"))
cat("Source Data Extended Data 3 (consolidated a-f) written:", paste(names(SHEETS), collapse = ", "), "\n")
