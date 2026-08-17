# =============================================================================
# Source Data — Figure 6 (secreted-exerkine cis-pQTL <-> GWAS disease network).
# Rebuilds the coloc network from the committed analysis inputs (res.olink.linear +
# the cis-pQTL coloc results + HPA secretome) via the shared network builder, then
# exports the panel tables. Independent of the figure-rendering script.
#
# Output: source_data/SourceData_Figure6.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(stringr); library(readr)
  library(readxl); library(tibble); library(igraph)
  library(clusterProfiler); library(org.Hs.eg.db); library(writexl) })
source(here("R/figures/helpers_coloc_network.R"))   # cn_prepare_coloc_traits, cn_build_trait_network_full
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

load(here("data", "res.olink.linear.rda"))
coloc_res     <- readxl::read_xlsx(here("data-raw", "phewas_pqtl", "newest_coloc_results_resubmission.xlsx"))
all_traits    <- read.csv(here("data-raw", "phewas_pqtl", "final_pqtl_gwas_with_corrected_parents.csv"))
trait_colors  <- c(
  "Digestive system disorder" = "#B7704C", "Cardiovascular disease" = "#b33232",
  "Metabolic disorder" = "#fdb462", "Immune system disorder" = "#ffed6f",
  "Neurological disorder" = "#ffffb3", "Lipid or lipoprotein measurement" = "#b3de69",
  "Inflammatory measurement" = "#ccebc5", "Hematological measurement" = "#8dd3c7",
  "Body measurement" = "#66ccff", "Cardiovascular measurement" = "#80b1d3",
  "Other measurement" = "#006699", "Response to drug" = "#fccde5",
  "Biological process" = "#bebada", "Cancer" = "#bc80bd",
  "Other trait" = "#fb8072", "Other disease" = "#ff3399")
cluster_color_palette <- c("0" = "#B0B0B0", "1" = "#1F78B4", "2" = "#33A02C", "3" = "#E31A1C",
                           "4" = "#FDBF6F", "5" = "#FF7F00", "6" = "#6A3D9A", "7" = "#CAB2D6")
hpa_clean <- readr::read_tsv(here("data-raw", "hpa_24.tsv"), show_col_types = FALSE) %>%
  dplyr::rename(gene = Gene, secretome_function = `Secretome function`) %>%
  dplyr::mutate(gene = str_trim(as.character(gene)),
                secretome_function = dplyr::na_if(str_squish(as.character(secretome_function)), ""),
                secretome_function = dplyr::if_else(str_to_lower(secretome_function) %in%
                  c("na", "n/a", "not available", "-", "none"), NA_character_, secretome_function))

disease_terms <- c("Cardiovascular disease", "Metabolic disorder", "Immune system disorder",
                   "Cancer", "Digestive system disorder", "Neurological disorder")
disease_combined <- cn_prepare_coloc_traits(all_traits, coloc_res, res.olink.linear, hpa_clean) %>%
  dplyr::filter(parent_term %in% disease_terms)

# rebuild the network (PDFs go to a scratch prefix so figure outputs are untouched)
res6 <- cn_build_trait_network_full(
  df = disease_combined, olink_df = res.olink.linear,
  out_prefix = file.path(tempdir(), "figure6_sourcedata_net"),
  title_text = "Secreted Exerkine cis pQTL <-> Open GWAS Disease Network",
  trait_palette = trait_colors, color_col = "parent_term", use_secreted = TRUE,
  min_mod_plot = 10L, alpha_adj = 0.05, run_fisher = TRUE, fisher_col = "parent_term",
  top_n_proteins = 10, top_n_traits = 5, cluster_color_palette = cluster_color_palette)

SHEETS[["b_network_nodes"]] <- igraph::as_data_frame(res6$g, "vertices")
SHEETS[["b_network_edges"]] <- igraph::as_data_frame(res6$g, "edges")
SHEETS[["c_top_traits"]]    <- res6$trait_cluster_comp
SHEETS[["d_top_proteins"]]  <- res6$prot_trait_comp
if (!is.null(res6$fisher_df)) SHEETS[["e_module_trait_enrichment"]] <- res6$fisher_df

# ---- f: unique GO:BP enriched terms per module (committed ORA workbook) ---------
ora_path <- here("multiomics_II", "figures", "network_diseases_secreted_testing_Combined_ORA.xlsx")
if (file.exists(ora_path)) {
  SHEETS[["f_unique_go"]] <- readxl::read_xlsx(ora_path) %>%
    dplyr::group_by(Module) %>% dplyr::arrange(p.adjust, desc(Count)) %>%
    dplyr::slice_head(n = 2) %>% dplyr::ungroup() %>%
    dplyr::mutate(log10p = -log10(p.adjust)) %>%
    dplyr::select(Module, Description, Count, p.adjust, log10p) %>% as.data.frame()
}

writexl::write_xlsx(SHEETS, here("source_data", "SourceData_Figure6.xlsx"))
cat("Source Data Figure 6 written:", paste(names(SHEETS), collapse = ", "), "\n")
