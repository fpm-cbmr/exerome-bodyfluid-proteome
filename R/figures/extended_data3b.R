# =============================================================================
# Extended Data Fig. 3 d-f (related to main Fig. 6) — Human genetic links between exercise-regulated exerkines and disease risk, and exerkines influenced by exercise modality and physiological traits
#
# Colocalization of the exercise-MODE-regulated secreted plasma proteins
# (replication cohort; group main-effect FDR < 0.05) with non-disease ("other")
# GWAS traits; trait nodes coloured by organ-system `category`. Built with
# R/figures/helpers_coloc_network.R (set.seed(123)).
#
# Panels:
#   d  network   e  top proteins (by connected traits)
#   f  top traits (by connected-protein exercise mode: A/H/S, or "Multiple")
#
# NOTE (panel f): each protein's exercise mode is derived from combined.pilot.res.linear
#   as the group(s) A/H/S in which its per-mode time model is significant, collapsed to
#   "Multiple" when >1. This is an interpretation — confirm the intended rule if panel f
#   must match the published version exactly (adjust mode_map below; nothing else changes).
#
# Inputs:
#   data/res.olink.linear.rda                                <- 04_plasma_lmm_clusters.R
#   data/res.validation.linear.pilot.rda,
#     data/combined.pilot.res.linear.rda                     <- 13_package_replication_period.R
#   data/validation_effect_size_clusters.rda                 <- 13b_replication_figure_objects.R
#   data/validation_prot.label.rda                           <- 10b_replication_olink_processing.R
#   data-raw/phewas_pqtl/*.{xlsx,csv}, data-raw/hpa_24.tsv    (coloc results + trait parents + HPA)
#
# Output: figures/extended_data_3/*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr); library(readr); library(readxl)
  library(igraph); library(ggplot2); library(ggpubr); library(writexl); library(ggsci)
})
source(here("R/figure_defaults.R"))
source(here("R/figures/helpers_coloc_network.R"))
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate; rename <- dplyr::rename
FIG <- here("figures", "extended_data_3"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(here("source_data"), showWarnings = FALSE)
OUT <- file.path(FIG, "exercise_mode_network_other_traits"); SHEETS <- list()

# ---- inputs ------------------------------------------------------------------
load(here("data", "res.olink.linear.rda"))
load(here("data", "validation_prot.label.rda")); prot.label <- validation_prot.label
load(here("data", "res.validation.linear.pilot.rda"))
res_grouxtime0 <- res.validation.linear.pilot
grp_clusters <- get(load(here("data", "validation_effect_size_clusters.rda")))
combined.pilot.res.linear <- get(load(here("data", "combined.pilot.res.linear.rda")))
coloc_res     <- read_xlsx(here("data-raw", "phewas_pqtl", "newest_coloc_results_resubmission.xlsx"))
parent_traits <- read.csv(here("data-raw", "phewas_pqtl", "final_pqtl_gwas_with_corrected_parents.csv"))
hpa_clean <- read_tsv(here("data-raw", "hpa_24.tsv"), show_col_types = FALSE) %>%
  rename(gene = Gene, secretome_function = `Secretome function`) %>%
  mutate(gene = str_trim(as.character(gene)), secretome_function = str_squish(as.character(secretome_function)),
         secretome_function = na_if(secretome_function, ""),
         secretome_function = if_else(str_to_lower(secretome_function) %in%
                                        c("na", "n/a", "not available", "-", "none"),
                                      NA_character_, secretome_function))

# ---- res_grouxtime: mode-significant proteins + replication cluster ----------
res_grouxtime <- res_grouxtime0 %>%
  rename(Assay = outcome) %>% filter(pval.aov.group.adj < 0.05) %>%
  inner_join(grp_clusters, by = "Assay") %>% rename(cluster = Cluster)

# prepare_coloc_traits_validation: harmonize traits, proteins, and HPA labels
prepare_coloc_traits_validation <- function(trait_df, coloc_res, res_valid, hpa_clean) {
  trait_map <- trait_df %>%
    rename(reported_traits = dplyr::any_of(c("reported_traits", "reported_trait"))) %>%
    select(reported_traits, parent_term, category) %>%
    mutate(reported_traits = str_trim(as.character(reported_traits)))
  coloc_df <- coloc_res %>%
    mutate(left_trait = str_trim(as.character(left_trait)),
           right_geneSymbol = str_to_upper(str_trim(as.character(right_geneSymbol)))) %>%
    inner_join(trait_map, by = c("left_trait" = "reported_traits"))
  olink_cluster <- res_valid %>%
    mutate(Assay = str_to_upper(str_trim(as.character(Assay)))) %>%
    select(Assay, cluster, beta.exposure0) %>% distinct(Assay, .keep_all = TRUE)
  coloc_df %>% inner_join(olink_cluster, by = c("right_geneSymbol" = "Assay")) %>%
    left_join(hpa_clean %>% select(gene, secretome_function) %>%
                mutate(gene = str_to_upper(str_trim(as.character(gene)))) %>% distinct(),
              by = c("right_geneSymbol" = "gene"))
}

other_traits   <- parent_traits %>% filter(grepl("^Other", parent_term, ignore.case = TRUE))
other_combined <- prepare_coloc_traits_validation(other_traits, coloc_res, res_grouxtime, hpa_clean)
categories <- sort(unique(other_combined$category))
category_palette <- setNames(ggsci::pal_d3("category20")(length(categories)), categories)

# ---- panels d & e via the shared builder ------------------------------------
res3b <- cn_build_trait_network_full(
  df = other_combined, olink_df = prot.label, out_prefix = OUT,
  title_text = "Secreted Exercise Mode Exerkines <-> Other Traits Network",
  trait_palette = category_palette, color_col = "category", use_secreted = TRUE,
  min_mod_plot = 10L, alpha_adj = 0.05, run_fisher = TRUE, fisher_col = "category",
  top_n_proteins = 10, top_n_traits = 5, emit_cluster_bar = FALSE)  # panel f is by mode, built below
SHEETS[["d_network_edges"]] <- igraph::as_data_frame(res3b$g, "edges")
SHEETS[["e_top_proteins"]]  <- res3b$prot_trait_comp

# ---- panel f: top traits coloured by connected-protein EXERCISE MODE ---------
# mode palette plus a "Multiple" category for mixed-mode connections
mode_pal <- c(A = "#004488", H = "#BB5566", S = "#DDAA33", Multiple = "grey70")
mode_map <- combined.pilot.res.linear %>%
  mutate(sig = (pval.exposure0.adj < 0.05) | (pval.exposure0.5.adj < 0.05) | (pval.exposure24.adj < 0.05)) %>%
  filter(sig) %>% mutate(Assay = str_to_upper(str_trim(as.character(outcome))),
                         group = recode(group, U = "A")) %>%
  group_by(Assay) %>% summarise(nmode = n_distinct(group), mode1 = dplyr::first(group), .groups = "drop") %>%
  mutate(Exercise_Mode = ifelse(nmode > 1, "Multiple", mode1))
edges_group <- res3b$edges_group
top_traits <- res3b$top_traits
f_df <- edges_group %>% filter(Trait %in% top_traits$name) %>%
  left_join(mode_map %>% select(Assay, Exercise_Mode), by = c("Protein" = "Assay")) %>%
  mutate(Exercise_Mode = ifelse(is.na(Exercise_Mode), "Multiple", Exercise_Mode)) %>%
  group_by(Trait, Exercise_Mode) %>% summarise(n = dplyr::n(), .groups = "drop") %>%
  group_by(Trait) %>% mutate(total_n = sum(n)) %>% ungroup() %>%
  mutate(Trait = str_wrap(Trait, 25))
pf <- ggplot(f_df, aes(reorder(Trait, total_n), n, fill = Exercise_Mode)) +
  geom_col(position = "stack") + coord_flip() +
  scale_fill_manual(values = mode_pal, na.value = "grey70") + theme_pubr(base_size = 6) +
  labs(title = "Top Traits: composition of connected protein clusters",
       x = "Trait", y = "Connections", fill = "Exercise Mode") +
  theme(text = element_text(size = 6), legend.position = "right",
        legend.key.size = unit(2, "mm"), panel.grid = element_blank()) + theme_strokes
ggsave(paste0(OUT, "_TopTraits_ExerciseMode.pdf"), pf, width = 7, height = 5, units = "cm", dpi = 600)
SHEETS[["f_top_traits_mode"]] <- f_df

# source data now built by R/source_data/_source_data.R
cat("Extended Data 3 (d-f) panels written:", paste(names(SHEETS), collapse = ", "), "\n")
