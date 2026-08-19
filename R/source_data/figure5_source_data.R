# =============================================================================
# Source Data — Figure 5 (replication cohort).
#
# Input:  data/res.olink.linear.rda, data/res.validation.linear.pilot.rda,
#         data/agree_df.rda, data/betas_sig.rda,
#         data/validation_effect_size_clusters.rda,
#         data/heatmap_timegroup_significant_cluster_enrichment.rda,
#         results/ora/concordant_ORA_all.csv
# Output: source_data/SourceData_Figure5.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(stringr); library(tibble); library(readr); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()
lev <- c("0h", "0.5h", "24h")

# ---- b: Venn set sizes (discovery time / replication time / replication mode) -
load(here("data", "res.olink.linear.rda"))
load(here("data", "res.validation.linear.pilot.rda"))
exerome_proteins    <- res.olink.linear %>% dplyr::filter(fdr.aov < 0.05) %>% dplyr::pull(Assay) %>% unique()
validation_proteins <- res.validation.linear.pilot %>% dplyr::filter(pval.aov.t.factor.adj < 0.05) %>% dplyr::pull(outcome) %>% unique()
validation_group    <- res.validation.linear.pilot %>%
  dplyr::filter(pval.aov.group.adj < 0.05 & pval.aov.t.factor.adj > 0.05) %>% dplyr::pull(outcome) %>% unique()
SHEETS[["b_venn_counts"]] <- data.frame(
  Set = c("Discovery_Time", "Replication_Time", "Replication_Exercise_Mode"),
  N = c(length(exerome_proteins), length(validation_proteins), length(validation_group)))

# ---- c: proportion concordant per time ---------------------------------------
load(here("data", "agree_df.rda"))
SHEETS[["c_concordance"]] <- as.data.frame(agree_df)

# ---- d: beta(discovery) vs beta(replication) scatter -------------------------
load(here("data", "betas_sig.rda"))
SHEETS[["d_beta_scatter"]] <- betas_sig %>%
  dplyr::select(Assay, Time_Point, Beta_Exerome, Beta_Validation, sig_cat) %>% as.data.frame()

# ---- e: selected enrichment terms for concordant exerkines -------------------
ora_tidy <- readr::read_csv(here("results", "ora", "concordant_ORA_all.csv"), show_col_types = FALSE)
focus <- tibble::tribble(
  ~label,                                    ~pattern,
  "Muscle contraction",                      "(^|\\b)muscle\\s+contraction(\\b|$)",
  "CS/DS degradation",                       "cs/ds\\s+degradation|chondroitin\\s+sulfate/.?dermatan\\s+sulfate\\s+degradation",
  "platelet activation",                     "(^|\\b)platelet\\s+activation(\\b|$)",
  "translation",                             "(^|\\b)translation(\\b|$)",
  "Signal Transduction",                     "(^|\\b)signal\\s+transduction(\\b|$)",
  "Metabolism of lipids",                    "metabolism\\s+of\\s+lipids(\\s+and\\s+lipoproteins)?",
  "integrin-mediated signaling pathway",     "integrin[-\\s]mediated\\s+signaling\\s+pathway",
  "vesicle-mediated transport",              "vesicle[-\\s]mediated\\s+transport",
  "translational initiation",                "translational\\s+initiation",
  "TNFs bind their physiological receptors", "tnfs\\s+bind\\s+their\\s+physiological\\s+receptors",
  "Extracellular matrix organization",       "extracellular\\s+matrix\\s+organization",
  "cytokine production",                     "(^|\\b)cytokine\\s+production(\\b|$)",
  "neutrophil activation",                   "neutrophil\\s+activation",
  "bone resorption",                         "bone\\s+resorption",
  "glial cell activation",                   "glial\\s+cell\\s+activation",
  "lipid transport",                         "lipid\\s+transport",
  "amyloid-beta clearance",                  "amyloid[-\\s]?beta\\s+clearance",
  "GPCR ligand binding",                     "gpcr\\s+ligand\\s+binding",
  "synaptic signaling",                      "synaptic\\s+signaling",
  "neuron development",                      "neuron\\s+development",
  "gastric motility",                        "gastric\\s+motility")
norm <- function(x) stringr::str_to_lower(stringr::str_squish(x))
ora_clean <- ora_tidy %>% dplyr::filter(!is.na(term_name)) %>%
  dplyr::mutate(term_name_norm = norm(term_name), Time_Point = factor(Time_Point, levels = lev),
                Direction = dplyr::recode(Direction, pos = "Positive", neg = "Negative"),
                dir_short = dplyr::if_else(Direction == "Positive", "(+)", "(-)"))
SHEETS[["e_concordant_terms"]] <- tidyr::crossing(ora_clean, focus) %>%
  dplyr::filter(stringr::str_detect(term_name_norm, stringr::regex(pattern, ignore_case = TRUE))) %>%
  dplyr::mutate(label = factor(label, levels = focus$label)) %>%
  dplyr::arrange(Time_Point, Direction, label, p_value) %>%
  dplyr::distinct(Time_Point, Direction, label, .keep_all = TRUE) %>%
  dplyr::mutate(gene_ratio = intersection_size / query_size, color_score = -log10(p_value)) %>%
  dplyr::select(label, Time_Point, dir_short, gene_ratio, color_score, term_name) %>% as.data.frame()

# ---- f: exercise-mode heatmap clusters (10 k-means, seed 123 = analysis 13b) --
load(here("data", "validation_effect_size_clusters.rda"))
SHEETS[["f_mode_heatmap_clusters"]] <- as.data.frame(validation_effect_size_clusters)

# ---- g: selected ORA terms per replication cluster ---------------------------
load(here("data", "heatmap_timegroup_significant_cluster_enrichment.rda"))
target_terms <- c("regeneration", "positive regulation of cell development",
  "extracellular matrix organization", "positive regulation of osteoblast differentiation",
  "response to oxidative stress", "organelle disassembly", "leukocyte activation", "cytokine production",
  "somatodendritic compartment", "coated vesicle", "synapse organization",
  "insulin secretion (glucose stimulus)", "neutrophil mediated cytotoxicity", "hormone secretion")
SHEETS[["g_cluster_ora"]] <- heatmap_timegroup_significant_cluster_enrichment %>%
  dplyr::mutate(Description = ifelse(Description == "insulin secretion involved in cellular response to glucose stimulus",
                                     "insulin secretion (glucose stimulus)", Description)) %>%
  dplyr::filter(Description %in% target_terms) %>%
  dplyr::select(cluster, Description, Count, p.adjust) %>% as.data.frame()

writexl::write_xlsx(SHEETS, here("source_data", "SourceData_Figure5.xlsx"))
cat("Source Data Figure 5 written:", paste(names(SHEETS), collapse = ", "), "\n")
