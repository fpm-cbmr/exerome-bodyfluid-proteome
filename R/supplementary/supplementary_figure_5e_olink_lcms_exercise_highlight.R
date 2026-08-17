# ==============================================================================
# Supplementary Figure 5e - Plasma proteome response to acute exercise measured with LC-MS/MS based proteomics (related to Fig. 3)
# ------------------------------------------------------------------------------
# Of the 170 plasma proteins measured on both platforms, which are exercise-
# regulated in the Olink discovery model (time main effect, fdr.aov < 0.05), and
# for those, is the Olink-vs-MS correlation positive or negative?
#
# Output:
#   - figures/supplementary/S5e_olink_lcms_correlation.csv
#   - figures/supplementary/S5e_olink_lcms_correlation_highlight.pdf
# ==============================================================================

suppressMessages({library(here); library(dplyr); library(tidyr)
                  library(ggplot2); library(ggpubr); library(ggrepel)})

load(here("data/plasma_npx_data.rda"))   # LC-MS plasma
load(here("data/olink_npx.data.rda"))    # Olink (OlinkID cols)
load(here("data/res.olink.linear.rda"))

# rename Olink OIDs -> gene symbols
map <- res.olink.linear %>% select(OlinkID, Assay) %>% distinct()
oid <- colnames(olink_npx.data)
ren <- map$Assay[match(oid, map$OlinkID)]; ren[is.na(ren)] <- oid[is.na(ren)]
colnames(olink_npx.data) <- ren
olink_npx.data <- olink_npx.data[, !duplicated(colnames(olink_npx.data))]

common <- setdiff(intersect(colnames(plasma_npx_data), colnames(olink_npx.data)),
                  c("label", "time", "subject"))

pl <- plasma_npx_data %>% select(label, time, any_of(common)) %>%
  pivot_longer(-c(label, time), names_to = "protein", values_to = "ms") %>%
  mutate(time = as.character(time))
ol <- olink_npx.data %>% select(label, time, any_of(common)) %>%
  pivot_longer(-c(label, time), names_to = "protein", values_to = "olink") %>%
  mutate(time = as.character(time))
m <- inner_join(pl, ol, by = c("label", "time", "protein")) %>%
  filter(is.finite(ms), is.finite(olink))

# per-protein Olink-vs-MS Spearman r, p, and MS abundance
res <- m %>% group_by(protein) %>%
  filter(n() >= 4) %>%
  summarise(n = n(),
            r = cor(ms, olink, method = "spearman"),
            r_p = suppressWarnings(cor.test(ms, olink, method = "spearman")$p.value),
            ms_abundance = mean(ms, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(r_fdr = p.adjust(r_p, "BH"))

# Olink exercise-regulated status (time main effect, fdr.aov < 0.05)
ex <- res.olink.linear %>% filter(fdr.aov < 0.05) %>% pull(Assay) %>% unique()
res <- res %>%
  mutate(exercise_reg = protein %in% ex,
         category = case_when(
           exercise_reg & r >= 0 ~ "Exercise-reg, positive r",
           exercise_reg & r <  0 ~ "Exercise-reg, negative r",
           TRUE ~ "Not exercise-reg"),
         category = factor(category, levels = c("Exercise-reg, positive r",
                                                "Exercise-reg, negative r",
                                                "Not exercise-reg")))

write.csv(res %>% arrange(category, r),
          here("figures/supplementary/S5e_olink_lcms_correlation.csv"),
          row.names = FALSE)

n_ex   <- sum(res$exercise_reg)
n_pos  <- sum(res$exercise_reg & res$r >= 0)
n_neg  <- sum(res$exercise_reg & res$r < 0)
n_negs <- sum(res$exercise_reg & res$r < 0 & res$r_fdr < 0.05)
n_poss <- sum(res$exercise_reg & res$r > 0 & res$r_fdr < 0.05)
cat(sprintf("Common proteins: %d | Olink exercise-regulated: %d\n", nrow(res), n_ex))
cat(sprintf("  positive r: %d (%d at FDR<0.05) | negative r: %d (%d at FDR<0.05)\n",
            n_pos, n_poss, n_neg, n_negs))
cat(sprintf("  median r among exercise-regulated: %.2f\n",
            median(res$r[res$exercise_reg])))

# labels: all negatively-correlated exercise-regulated proteins + all strong
# positives (exercise-regulated with r >= 0.5, i.e. above the dashed line)
lab <- bind_rows(
  res %>% filter(exercise_reg, r < 0) %>% arrange(r),
  res %>% filter(exercise_reg, r >= 0.5) %>% arrange(desc(r)))
cat(sprintf("labelling %d proteins (%d positive r>=0.5, %d negative)\n",
            nrow(lab), sum(lab$r >= 0.5), sum(lab$r < 0)))

cols <- c("Exercise-reg, positive r" = "#B2182B",
          "Exercise-reg, negative r" = "#2166AC",
          "Not exercise-reg" = "grey80")

# --- unit conversions: true points ------------------------------------------
PT   <- 72.27 / 25.4        # ggplot's mm<->pt factor (.pt)
lw   <- 0.5 / PT            # 0.5 pt stroke lines (ggplot linewidth is in mm)
fs   <- 6   / PT            # 6 pt for geom_text/label (size is in mm)
ff   <- "Arial"

p <- ggplot(res, aes(x = ms_abundance, y = r)) +
  geom_hline(yintercept = 0, linewidth = lw, color = "grey50") +
  geom_hline(yintercept = 0.5, linewidth = lw, linetype = "dashed", color = "grey70") +
  geom_point(data = filter(res, !exercise_reg), color = "grey80", size = 0.7,
             alpha = 0.6, stroke = lw) +
  geom_point(data = filter(res, exercise_reg), aes(color = category), size = 1.2,
             alpha = 0.9, stroke = lw) +
  ggrepel::geom_text_repel(data = lab, aes(label = protein, color = category), size = fs,
                           family = ff, min.segment.length = 0, segment.size = lw,
                           max.overlaps = Inf, box.padding = 0.15, point.padding = 0.1,
                           force = 2, show.legend = FALSE) +
  scale_color_manual(values = cols, name = NULL, breaks = names(cols)[1:2]) +
  labs(x = "Mean LC-MS abundance (log2)",
       y = "Olink vs LC-MS correlation (Spearman r)",
       title = "Cross-platform agreement of Olink exercise-regulated proteins",
       subtitle = sprintf("%d/%d shared proteins exercise-regulated: %d positive r, %d negative r",
                          n_ex, nrow(res), n_pos, n_neg)) +
  theme_pubr(base_size = 6, base_family = ff, legend = "right") +
  theme(text = element_text(size = 6, family = ff),
        plot.title = element_text(size = 6, face = "bold", family = ff),
        plot.subtitle = element_text(size = 6, family = ff),
        axis.title = element_text(size = 6, family = ff),
        axis.text = element_text(size = 6, family = ff),
        legend.text = element_text(size = 6, family = ff),
        legend.title = element_text(size = 6, family = ff),
        legend.key.size = unit(2.5, "mm"),
        axis.line = element_line(linewidth = lw),
        axis.ticks = element_line(linewidth = lw))

ggsave(here("figures/supplementary/S5e_olink_lcms_correlation_highlight.pdf"),
       p, width = 11, height = 5, units = "cm", dpi = 600, device = cairo_pdf)

source(here("R", "supplementary", "helpers_supp_plots.R"))
# source data now built by R/source_data/ (per-figure script)   # merges with a-d
cat("\nWrote figures/supplementary/S5e_olink_lcms_correlation.csv and",
    "figures/supplementary/S5e_olink_lcms_correlation_highlight.pdf\n")
