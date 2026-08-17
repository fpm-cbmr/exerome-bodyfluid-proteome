# =============================================================================
# Figure 1 — Acute exercise triggers compartmentalized temporal dynamics
# across body fluids.
#
# Panels:
#   a  study-design schematic (BioRender; no source data)
#   b  plasma metabolic markers over time
#   c  MSD cytokine responses over time
#   d  no. significantly regulated proteins per fluid / time point
#   e  UpSet of proteins shared across >=2 fluids, by directional concordance
#
# Inputs (produced by the analysis pipeline):
#   data/res.msd.linear.rda, data/msd.exerome.dat.rda,
#     data/metabolic.cytokine.prot.label.rda   <- 01_metabolic_cytokine_lmm.R   (b,c)
#   data/res.saliva.linear.rda                  <- 02_saliva_lmm_clusters.R      (d,e)
#   data/res.urine.linear.rda                   <- 03_urine_lmm_clusters.R       (d,e)
#   data/res.olink.linear.rda                   <- 04_plasma_lmm_clusters.R      (d,e)
#
# Output: figures/figure_1/*.pdf
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggpubr); library(patchwork); library(writexl)
})
source(here("R/package_loading.R"))
source(here("R/functions_loading.R"))
source(here("R/figure_defaults.R"))   # 0.5 pt strokes throughout
# pin dplyr verbs (package_loading/functions_loading attach packages that mask them)
count <- dplyr::count; filter <- dplyr::filter; select <- dplyr::select
mutate <- dplyr::mutate; summarise <- dplyr::summarise
rename <- dplyr::rename; arrange <- dplyr::arrange
dir.create(here("figures", "figure_1"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("source_data"),         showWarnings = FALSE)

# ---- inputs (saved analysis outputs) ----------------------------------------
load(here("data", "res.msd.linear.rda"))                 # from analysis/01
load(here("data", "msd.exerome.dat.rda"))
load(here("data", "metabolic.cytokine.prot.label.rda"))
load(here("data", "res.olink.linear.rda"))               # from analysis/04
load(here("data", "res.saliva.linear.rda"))              # from analysis/02
load(here("data", "res.urine.linear.rda"))               # from analysis/03
prot.label <- metabolic.cytokine.prot.label

# ---- panels b, c: metabolic markers & cytokines -----------------------------
cl.label.exerome <- tibble(cluster = 1:2,
                        label = c("Metabolic markers", "MDS cytokine panel"),
                        cluster.paper = c("", ""))
exerome.dat <- msd.exerome.dat %>% dplyr::rename(time_label = time)   # metabolic_markers() expects time_label

pb <- metabolic_markers(res.msd.linear, exerome.dat, prot.label, cl.label.exerome, 1) +
  theme(text = element_text(size = 6), legend.position = c(0.8, 1)) +
  coord_cartesian(ylim = c(-1, 3)) + theme_strokes
ggsave(here("figures", "figure_1", "b_metabolic_markers.pdf"), pb,
       width = 5.5, height = 5, units = "cm", dpi = 600)

pc <- metabolic_markers(res.msd.linear, exerome.dat, prot.label, cl.label.exerome, 2) +
  theme(text = element_text(size = 6)) + coord_cartesian(ylim = c(-0.5, 2)) + theme_strokes
ggsave(here("figures", "figure_1", "c_mds_cytokines.pdf"), pc,
       width = 5.5, height = 5, units = "cm", dpi = 600)

# ---- panel d: significant-protein counts per fluid/timepoint/direction -------
count_sig <- function(df, fluid_name) {
  pv <- grep("^pval\\.exposure", names(df), value = TRUE)
  for (p in pv) df[[paste0(p, "_fdr")]] <- p.adjust(df[[p]], "fdr")
  fdr <- df %>% select(Assay, matches("pval\\.exposure.*_fdr")) %>%
    pivot_longer(-Assay, names_to = "time_point", names_pattern = "pval\\.exposure(.*)_fdr", values_to = "pval_fdr")
  bet <- df %>% select(Assay, matches("^beta\\.exposure")) %>%
    pivot_longer(-Assay, names_to = "time_point", names_pattern = "beta\\.exposure(.*)", values_to = "beta")
  fdr %>% inner_join(bet, by = c("Assay", "time_point"), relationship = "many-to-many") %>%
    filter(pval_fdr < 0.05) %>%
    mutate(direction = ifelse(beta > 0, "increasing", "decreasing"), fluid = fluid_name) %>%
    group_by(fluid, time_point, direction) %>% summarise(count = n(), .groups = "drop")
}
panel_d <- bind_rows(count_sig(res.saliva.linear, "Saliva"),
                     count_sig(res.urine.linear,  "Urine"),
                     count_sig(res.olink.linear,  "Plasma")) %>%
  mutate(time_point = factor(time_point, levels = c("0", "0.5", "1", "3", "24")))
fluid_colors <- c("Plasma" = "#4B0082", "Saliva" = "#FFA07A", "Urine" = "#FFD700")
time_color_palette <- c("0" = "#56B4E9", "0.5" = "#009E73", "1" = "#F0E442", "3" = "#0072B2", "24" = "#CC79A7")
total_counts <- panel_d %>% group_by(time_point) %>% summarise(total = sum(count), .groups = "drop")

pd <- ggplot(panel_d, aes(x = time_point, y = count, fill = fluid)) +
  geom_bar(data = panel_d %>% filter(direction == "increasing"), aes(y = count),
           stat = "identity", position = position_dodge(width = 0.8), width = 0.75) +
  geom_bar(data = panel_d %>% filter(direction == "decreasing"), aes(y = -count),
           stat = "identity", position = position_dodge(width = 0.8), width = 0.75) +
  geom_hline(yintercept = 0, color = "black") +
  geom_text(data = total_counts,
            aes(x = time_point, y = max(panel_d$count) * 1.1, label = total, color = time_point),
            inherit.aes = FALSE, size = 2.5) +
  scale_fill_manual(values = fluid_colors) +
  scale_color_manual(values = time_color_palette) +
  scale_y_continuous(labels = abs) +
  theme_pubr(base_size = 6) +
  labs(x = "Time Point", y = "Counts of decreased(-)/increased(+) proteins",
       fill = "Fluid", color = "Time point") +
  theme(axis.text = element_text(size = 6), axis.title = element_text(size = 6),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        legend.key.size = unit(2, "mm"), legend.position = "right") + theme_strokes
ggsave(here("figures", "figure_1", "d_significant_proteins.pdf"), pd,
       width = 11, height = 5, units = "cm", dpi = 600)

# ---- panel e: UpSet of proteins shared across >=2 fluids, by direction --------
fluid_direction <- function(res, fluid_name) {
  res <- as.data.frame(res)
  beta_cols <- grep("^beta\\.exposure", colnames(res), value = TRUE)
  sig <- res[res$fdr.aov < 0.05 & !is.na(res$fdr.aov), ]
  peak <- apply(as.matrix(sig[, beta_cols, drop = FALSE]), 1, function(b) {
    b <- b[!is.na(b)]; if (length(b) == 0) NA_real_ else b[which.max(abs(b))]
  })
  tibble(Assay = sig$Assay, fluid = fluid_name,
         direction = ifelse(peak > 0, "up", "down")) %>%
    filter(!is.na(direction)) %>% distinct(Assay, .keep_all = TRUE)
}
dir_tab <- bind_rows(fluid_direction(res.olink.linear, "Plasma"),
                     fluid_direction(res.saliva.linear, "Saliva"),
                     fluid_direction(res.urine.linear, "Urine"))
prot <- dir_tab %>% group_by(Assay) %>%
  summarise(fluids = paste(sort(fluid), collapse = "+"), n_fluids = n(),
            n_up = sum(direction == "up"), n_down = sum(direction == "down"), .groups = "drop") %>%
  filter(n_fluids >= 2) %>%
  mutate(Direction = case_when(n_down == 0 ~ "Concordant up",
                               n_up == 0   ~ "Concordant down",
                               TRUE        ~ "Discordant"))
dir_cols <- c("Concordant up" = "#B2182B", "Concordant down" = "#2166AC", "Discordant" = "#999999")
int_levels <- prot %>% count(fluids) %>% arrange(desc(n)) %>% pull(fluids)
prot <- prot %>% mutate(fluids = factor(fluids, levels = int_levels))
bar_df <- prot %>% count(fluids, Direction) %>% mutate(Direction = factor(Direction, levels = names(dir_cols)))
tot_df <- prot %>% count(fluids)

p_top <- ggplot(bar_df, aes(fluids, n, fill = Direction)) +
  geom_col(width = 0.7, color = "white", linewidth = STROKE) +
  geom_text(data = tot_df, aes(fluids, n, label = n), inherit.aes = FALSE, vjust = -0.3, size = 1.9) +
  scale_fill_manual(values = dir_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) + labs(y = "Proteins") +
  theme_pubr(base_size = 6, legend = "none") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.text.y = element_text(size = 5),
        legend.position = "inside", legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1), legend.text = element_text(size = 5),
        legend.key.size = unit(2, "mm"), plot.margin = margin(1, 1, 0, 1)) + theme_strokes
fluid_rows <- c("Plasma", "Saliva", "Urine")
mat <- expand_grid(fluids = factor(int_levels, levels = int_levels),
                   fluid = factor(fluid_rows, levels = rev(fluid_rows))) %>%
  rowwise() %>% mutate(member = grepl(as.character(fluid), as.character(fluids), fixed = TRUE)) %>% ungroup()
seg <- mat %>% filter(member) %>% group_by(fluids) %>%
  summarise(ymin = min(as.integer(fluid)), ymax = max(as.integer(fluid)), .groups = "drop")
p_bot <- ggplot(mat, aes(fluids, fluid)) +
  geom_point(aes(color = member), size = 1.4) +
  geom_segment(data = seg, aes(x = fluids, xend = fluids, y = ymin, yend = ymax), inherit.aes = FALSE,
               linewidth = STROKE, color = "black") +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey85"), guide = "none") +
  labs(x = NULL, y = NULL) + theme_pubr(base_size = 6) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank(), plot.margin = margin(0, 1, 1, 1)) + theme_strokes
pe <- p_top / p_bot + plot_layout(heights = c(0.72, 0.28))
ggsave(here("figures", "figure_1", "e_upset_shared_directionality.pdf"), pe,
       width = 5.5, height = 5, units = "cm", dpi = 600, device = cairo_pdf)

# Source Data for this figure is built separately by
# R/source_data/figure1_source_data.R (writes source_data/SourceData_Figure1.xlsx).

cat(sprintf("Figure 1 done. Panel e: %d proteins in >=2 fluids (%s).\n",
            nrow(prot), paste(levels(prot$fluids), collapse = ", ")))
