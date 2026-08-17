# =============================================================================
# Supplementary Figure 6a — Supplementary Figure 6b — Detection of canonical exerkines and orthogonal validation of cytokines in the plasma Olink platform (related to Fig. 3)
#
# Quality-control check that well-established exercise-regulated proteins
# are detected on the Olink platform with their expected temporal
# dynamics. Combines a curated exerkine list (data-raw/exerkine_list.xlsx) with
# validated exerkines, and plots each protein's z-score at every
# post-exercise timepoint (* marks BH-FDR < 0.05 within a timepoint).
#
# Input:  data-raw/exerkine_list.xlsx, data/res.olink.linear.rda
# Output: figures/supplementary/S6a_exerkine_zscore_panel.pdf
#         figures/supplementary/S6a_exerkine_zscore_table.csv
# =============================================================================
suppressMessages({library(here); library(readxl); library(dplyr); library(tidyr)
                  library(stringr); library(ggplot2); library(ggpubr)})
OUT <- here("figures", "supplementary"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- formatting constants (true points) -------------------------------------
PT <- 72.27 / 25.4        # ggplot mm<->pt factor (.pt)
lw <- 0.5 / PT            # 0.5 pt strokes (linewidth is in mm)
fs <- 6   / PT            # 6 pt for geom_text (size is in mm)
ff <- "Arial"

# ---- exerkine lists: curated user list + literature-based canonical set -----
alias <- c(IL8 = "CXCL8", VEGF = "VEGFA")            # list name -> Olink Assay name
to_olink <- function(g) ifelse(g %in% names(alias), alias[g], g)

user_genes <- unique(toupper(trimws(read_excel(here("data-raw", "exerkine_list.xlsx"))$exerkine)))
user_olink <- unname(to_olink(user_genes))

mine_olink <- c("IL6","IL10","IL15","CXCL8","LIF","OSM","FGF21","GDF15","TNFSF11",
                "FABP4","GHRL","LEP","VEGFA","DCN","FST","ANGPTL4","BDNF","MSTN",
                "CCL2","SPP1")

all_olink <- union(user_olink, mine_olink)
gene_map <- tibble(olink = all_olink) %>%
  mutate(source = case_when(
    olink %in% user_olink & olink %in% mine_olink ~ "Both lists",
    olink %in% user_olink                          ~ "User list",
    TRUE                                           ~ "Curated (literature)"),
    list_name = olink)

# ---- discovery Olink model: z = beta / se per timepoint ---------------------
load(here("data", "res.olink.linear.rda"))
res <- as.data.frame(res.olink.linear)
res$.rid <- seq_len(nrow(res))
beta_cols <- grep("^beta\\.exposure", colnames(res), value = TRUE)
se_cols   <- grep("^se\\.exposure",   colnames(res), value = TRUE)
pv_cols   <- grep("^pval\\.exposure", colnames(res), value = TRUE)

long <- function(cols, val) res %>% select(.rid, Assay, all_of(cols)) %>%
  pivot_longer(-c(.rid, Assay), names_to = "tp",
               names_pattern = ".*exposure(.*)", values_to = val)

z_tab <- long(beta_cols, "beta") %>%
  inner_join(long(se_cols, "se") %>% select(.rid, tp, se), by = c(".rid", "tp")) %>%
  inner_join(long(pv_cols, "pval") %>% select(.rid, tp, pval), by = c(".rid", "tp")) %>%
  group_by(tp) %>% mutate(fdr = p.adjust(pval, "BH")) %>% ungroup() %>%
  mutate(z = beta / se)

panel <- z_tab %>%
  filter(Assay %in% gene_map$olink) %>%
  left_join(gene_map, by = c("Assay" = "olink")) %>%
  mutate(tp = factor(tp, levels = c("0", "0.5", "1", "3", "24")),
         sig = fdr < 0.05)

absent <- setdiff(gene_map$olink, unique(z_tab$Assay))
cat(sprintf("exerkines on Olink: %d / %d (absent: %s)\n",
            n_distinct(panel$Assay), nrow(gene_map), paste(absent, collapse = ", ")))
cat("counts by source:\n"); print(panel %>% distinct(list_name, source) %>% count(source))
write.csv(panel %>% select(source, list_name, Assay, tp, beta, se, z, pval, fdr, sig) %>%
            arrange(list_name, tp),
          file.path(OUT, "S6a_exerkine_zscore_table.csv"), row.names = FALSE)

# alphabetical order (A at top): reverse-sorted factor levels
panel <- panel %>% mutate(list_name = factor(list_name, levels = rev(sort(unique(list_name)))))

# ---- heatmap of z-scores ----------------------------------------------------
p <- ggplot(panel, aes(x = tp, y = list_name, fill = z)) +
  geom_tile(color = "white", linewidth = lw) +
  geom_text(data = filter(panel, sig), aes(label = "*"), size = fs, family = ff,
            vjust = 0.78, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "z-score") +
  labs(x = "Time post-exercise [h]", y = NULL,
       title = "Canonical exerkines (discovery plasma, Olink): z-score (beta/SE)") +
  theme_pubr(base_size = 6, base_family = ff, legend = "right") +
  theme(text = element_text(size = 6, family = ff),
        plot.title = element_text(size = 6, face = "bold", family = ff),
        axis.title = element_text(size = 6, family = ff),
        axis.text = element_text(size = 6, family = ff),
        legend.title = element_text(size = 6, family = ff),
        legend.text = element_text(size = 6, family = ff),
        legend.key.size = unit(3, "mm"),
        axis.line = element_line(linewidth = lw),
        axis.ticks = element_line(linewidth = lw))

ggsave(file.path(OUT, "S6a_exerkine_zscore_panel.pdf"),
       p, width = 11, height = 11, units = "cm", dpi = 600, device = cairo_pdf)
source(here("R", "supplementary", "helpers_supp_plots.R"))
# source data now built by R/source_data/ (per-figure script)   # merges with panel b (MSD)
cat("Wrote figures/supplementary/exerkine_zscore_panel.{pdf,csv}\n")
