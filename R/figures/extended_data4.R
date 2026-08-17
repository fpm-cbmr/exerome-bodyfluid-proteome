# =============================================================================
# Extended Data Fig. 4 (related to main Fig. 6) — Approved drug targets among the colocalized secreted exerkines.
#
# Panels (x = number of distinct approved DrugBank drug IDs per protein):
#   a  all secreted exercise-regulated exerkines with cis-pQTL trait coloc
#   b  secreted EXERCISE-MODE exerkines with cis-pQTL trait coloc
#
# Inputs:
#   data-raw/approved_drug_targets.csv                 (DrugBank approved targets)
#   data-raw/hpa_24.tsv                                (HPA secretome filter)
#   data/res.olink.linear.rda             <- 04_plasma_lmm_clusters.R          (a, discovery)
#   data/res.validation.linear.pilot.rda  <- 13_package_replication_period.R   (b, mode effect)
#
# Output: figures/extended_data_4/*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr); library(readr)
  library(ggplot2); library(ggpubr); library(forcats); library(writexl); library(snakecase)
})
source(here("R/figure_defaults.R"))
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate; rename <- dplyr::rename
FIG <- here("figures", "extended_data_4"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

# ---- approved drug targets (DrugBank) ---------------------------------------
adt <- read_csv(here("data-raw", "approved_drug_targets.csv"), show_col_types = FALSE)
colnames(adt) <- to_snake_case(colnames(adt))
adt <- adt %>% rename(drug_ids = drug_i_ds) %>% mutate(Assay = sub("_.*", "", uniprot_title))

# ---- secretome (HPA) --------------------------------------------------------
hpa_clean <- read_tsv(here("data-raw", "hpa_24.tsv"), show_col_types = FALSE) %>%
  rename(Assay = Gene, secretome_function = `Secretome function`) %>%
  mutate(Assay = str_trim(as.character(Assay)), secretome_function = str_squish(as.character(secretome_function)),
         secretome_function = na_if(secretome_function, ""),
         secretome_function = if_else(str_to_lower(secretome_function) %in%
                                        c("na", "n/a", "not available", "-", "none"),
                                      NA_character_, secretome_function))

# count distinct approved drugs per secreted, significant protein
drug_counts <- function(sig_assays) {
  adt %>% filter(Assay %in% sig_assays) %>%
    mutate(drug_ids = str_split(drug_ids, "[,;]\\s*")) %>% tidyr::unnest(drug_ids) %>%
    mutate(drug_ids = str_trim(drug_ids)) %>% filter(drug_ids != "") %>%
    distinct(Assay, drug_ids) %>%
    inner_join(hpa_clean, by = "Assay") %>% filter(!is.na(secretome_function)) %>%
    group_by(Assay) %>% summarise(n_drugs = n_distinct(drug_ids), .groups = "drop") %>%
    arrange(desc(n_drugs))
}
drug_bar <- function(df, ylab, out) {
  p <- ggplot(df, aes(fct_reorder(Assay, n_drugs), n_drugs)) +
    geom_col(fill = "#2171B5") + coord_flip() + theme_pubr(base_size = 6) +
    labs(title = NULL, x = ylab, y = "Number of approved drug targets") +
    theme(axis.text.y = element_text(size = 5), text = element_text(size = 6)) + theme_strokes
  ggsave(out, p, width = 8.25, height = 12, units = "cm", dpi = 600); p
}

# ---- panel a: all discovery secreted coloc exerkines ------------------------
load(here("data", "res.olink.linear.rda"))
a_assays <- res.olink.linear %>% filter(fdr.aov < 0.05) %>% pull(Assay) %>% unique()
a_df <- drug_counts(a_assays)
drug_bar(a_df, "Secreted exerkines that have cis pQTL and trait colocalization",
         file.path(FIG, "a_drug_targets_all.pdf"))
SHEETS[["a_drug_targets_all"]] <- a_df

# ---- panel b: secreted exercise-mode coloc exerkines ------------------------
load(here("data", "res.validation.linear.pilot.rda"))
b_assays <- res.validation.linear.pilot %>% filter(pval.aov.group.adj < 0.05) %>%
  pull(outcome) %>% unique()
b_df <- drug_counts(b_assays)
drug_bar(b_df, "Secreted exercise-mode exerkines that have cis pQTL and trait colocalization",
         file.path(FIG, "b_drug_targets_mode.pdf"))
SHEETS[["b_drug_targets_mode"]] <- b_df

# source data now built by R/source_data/_source_data.R
cat("Extended Data 4: panel a =", nrow(a_df), "proteins, panel b =", nrow(b_df), "proteins\n")
