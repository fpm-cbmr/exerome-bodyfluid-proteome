# =============================================================================
# Supplementary Figure 6b — Detection of canonical exerkines and orthogonal validation of cytokines in the plasma Olink platform (related to Fig. 3)
# The Discovery cohort has an independent MSD (Meso Scale Discovery) V-PLEX
# electrochemiluminescence immunoassay panel measured in the same participants
# and timepoints as Olink. Six cytokines are common to both platforms, giving an
# orthogonal validation of the Olink measurements with no new experiments:
#
#   MSD analyte  <->  Olink Assay
#   IL6          <->  IL6
#   IL8          <->  CXCL8
#   IL10         <->  IL10
#   TNFa         <->  TNF
#   IFNy         <->  IFNG
#   IL2          <->  IL2
#
# Both platforms use the same baseline mean-scaling; Spearman correlation is
# rank-based and so is unaffected by the per-analyte linear scaling.
#
# Input:  data/msd.exerome.dat.rda, data/olink.exerome.dat.rda, data/res.olink.linear.rda
# Output: figures/supplementary/S6b_msd_olink_validation.pdf
#         figures/supplementary/S6b_msd_olink_validation_stats.csv
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggpubr)
})
OUT <- here("figures", "supplementary"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

load(here("data", "msd.exerome.dat.rda"))     # MSD panel, baseline-scaled
load(here("data", "olink.exerome.dat.rda"))   # Olink, baseline-scaled (OID columns)
load(here("data", "res.olink.linear.rda"))    # OlinkID <-> Assay map

# harmonise the subject key (MSD stores "EX01", Olink stores "01") so the two
# platforms join on the same participant id
msd.exerome.dat$subject   <- sub("^EX", "", as.character(msd.exerome.dat$subject))
olink.exerome.dat$subject <- sub("^EX", "", as.character(olink.exerome.dat$subject))

# ---- shared cytokines and platform name mapping -----------------------------
shared <- tribble(
  ~label,    ~msd,    ~olink_assay,  ~plabel,
  "IL-6",    "IL6",   "IL6",         "'IL-6'",
  "IL-8",    "IL8",   "CXCL8",       "'IL-8'",
  "IL-10",   "IL10",  "IL10",        "'IL-10'",
  "TNF-a",   "TNFa",  "TNF",         "'TNF-'*alpha",
  "IFN-g",   "IFNy",  "IFNG",        "'IFN-'*gamma",
  "IL-2",    "IL2",   "IL2",         "'IL-2'"
)

oid_map <- res.olink.linear %>% select(OlinkID, Assay) %>% distinct() %>%
  filter(Assay %in% shared$olink_assay)

# ---- MSD long ---------------------------------------------------------------
msd_long <- msd.exerome.dat %>%
  select(subject, time, all_of(shared$msd)) %>%
  pivot_longer(-c(subject, time), names_to = "msd", values_to = "msd_value") %>%
  left_join(shared, by = "msd")

# ---- Olink long -------------------------------------------------------------
olink_sel <- olink.exerome.dat %>%
  select(subject, time = time_label, all_of(oid_map$OlinkID)) %>%
  pivot_longer(-c(subject, time), names_to = "OlinkID", values_to = "olink_value") %>%
  left_join(oid_map, by = "OlinkID") %>%
  left_join(shared %>% select(label, olink_assay), by = c("Assay" = "olink_assay"))

# ---- join same participant x timepoint x cytokine ---------------------------
paired <- inner_join(
  msd_long  %>% select(subject, time, label, msd_value),
  olink_sel %>% select(subject, time, label, olink_value),
  by = c("subject", "time", "label")
) %>% filter(is.finite(msd_value), is.finite(olink_value))

# ---- per-cytokine Spearman (all subject x timepoint samples) ----------------
stats <- paired %>%
  group_by(label) %>%
  summarise(
    n   = n(),
    rho = cor(msd_value, olink_value, method = "spearman", use = "complete.obs"),
    p   = suppressWarnings(cor.test(msd_value, olink_value, method = "spearman")$p.value),
    .groups = "drop"
  ) %>%
  arrange(desc(rho))

# overall pooled correlation across all cytokines/samples
overall <- with(paired, cor.test(msd_value, olink_value, method = "spearman"))
stats_all <- bind_rows(
  stats,
  tibble(label = "ALL (pooled)", n = nrow(paired),
         rho = unname(overall$estimate), p = overall$p.value)
)
write.csv(stats_all, file.path(OUT, "S6b_msd_olink_validation_stats.csv"), row.names = FALSE)

cat("\n===== MSD (immunoassay) vs Olink (PEA) cross-platform Spearman =====\n")
print(as.data.frame(stats_all), digits = 3)
cat(sprintf("\nMedian per-cytokine rho: %.2f | %d/%d cytokines rho>0.5 (p<0.05)\n",
            median(stats$rho), sum(stats$rho > 0.5 & stats$p < 0.05), nrow(stats)))

# ---- figure: faceted scatter with Spearman R + p ----------------------------
# plain text facet labels (ordered by rho) so every cytokine name uses the same
# font — no plotmath/Greek so fonts stay consistent across facets
paired <- paired %>% mutate(label = factor(as.character(label), levels = stats$label))

p_scatter <- ggplot(paired, aes(x = olink_value, y = msd_value)) +
  geom_point(size = 0.5, alpha = 0.6, color = "#4B0082") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, color = "#8B0000") +
  facet_wrap(~label, scales = "free", ncol = 3) +
  stat_cor(method = "spearman", size = 1.8, label.sep = "\n", cor.coef.name = "rho") +
  labs(x = "Olink NPX (baseline-scaled)",
       y = "MSD V-PLEX (baseline-scaled)",
       title = "Orthogonal validation: MSD immunoassay vs Olink") +
  theme_pubr(base_size = 6) +
  theme(strip.text = element_text(size = 6, face = "bold"),
        plot.title = element_text(size = 6, face = "bold"),
        axis.text  = element_text(size = 5))

ggsave(file.path(OUT, "S6b_msd_olink_validation.pdf"),
       p_scatter, width = 11, height = 7.5, units = "cm", dpi = 600, device = cairo_pdf)
source(here("R", "supplementary", "helpers_supp_plots.R"))
# source data now built by R/source_data/ (per-figure script)
cat("\nWrote figures/supplementary/S6b_msd_olink_validation.{pdf,csv}\n")
