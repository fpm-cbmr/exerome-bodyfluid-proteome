# =============================================================================
# Analysis 01 — Metabolic & cytokine panel: processing + mixed-effects models
# (Discovery cohort).
#
# Reads the raw MSD metabolic-marker / glucose / cytokine panel, reshapes it to a
# wide per-subject NPX table, builds the analyte annotation (metabolic vs cytokine
# cluster), baseline (t = -1) z-normalises, then fits a linear mixed model
# (time as fixed effect, subject random intercept, age + VO2max covariates)
# to each analyte.
#
# Input:  data-raw/clinical_data/glucose_cytokine_exercise_response.xlsx,
#         data-raw/clinical_data/covariates_msd.xlsx
# Output: data/metabolic.cytokine.prot.label.rda, data/msd.exerome.dat.rda,
#         data/res.msd.linear.rda   (consumed by R/figures/figure1.R, panels b,c)
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(readxl)
  library(data.table); library(doParallel)
})
source(here("R/package_loading.R"))
source(here("R/functions_loading.R"))

# ---- raw metabolic-marker / glucose / cytokine panel -> wide NPX table -------
data <- readxl::read_excel(here("data-raw", "clinical_data", "glucose_cytokine_exercise_response.xlsx"))
long_data <- data %>%
  tidyr::pivot_longer(cols = c(Pre, T00, T30, T60, T180, `24h`),
                      names_to = "time", values_to = "value") %>%
  dplyr::mutate(time = dplyr::case_when(
    time == "Pre" ~ -1, time == "T00" ~ 0,  time == "T30" ~ 0.5,
    time == "T60" ~ 1,  time == "T180" ~ 3, time == "24h" ~ 24),
    value = as.numeric(value))
metabolic.cytokine.npx.data <- long_data %>%
  tidyr::pivot_wider(names_from = analyte, values_from = value)
npx.data <- metabolic.cytokine.npx.data

# analyte annotation: cluster 1 = metabolic markers, cluster 2 = cytokines
metabolic.cytokine.prot.label <- long_data %>%
  dplyr::distinct(analyte) %>% dplyr::mutate(Assay = analyte) %>%
  dplyr::mutate(cluster = dplyr::case_when(
    Assay %in% c("CK_U.L", "c_peptide_pmol.L", "glucose_mmol.L", "insulin_mmol.L") ~ "1",
    TRUE ~ "2"))
save(metabolic.cytokine.prot.label, file = here("data", "metabolic.cytokine.prot.label.rda"))
prot.label <- metabolic.cytokine.prot.label

# ---- baseline (t = -1) z-normalisation ---------------------------------------
exerome.norm <- npx.data %>% dplyr::filter(time == -1) %>%
  tidyr::pivot_longer(cols = -c(time, subject), names_to = "Assay", values_to = "value") %>%
  dplyr::group_by(Assay) %>%
  dplyr::summarize(mean.value = mean(value, na.rm = TRUE), sd.value = sd(value, na.rm = TRUE), .groups = "drop")
exerome.dat <- npx.data
for (j in seq_len(nrow(exerome.norm)))
  exerome.dat <- exerome.dat %>%
    dplyr::mutate(!!exerome.norm$Assay[j] := (!!sym(exerome.norm$Assay[j]) - exerome.norm$mean.value[j]) / exerome.norm$sd.value[j])

# clinical covariates + factor time
clinical_data <- readxl::read_xlsx(here("data-raw", "clinical_data", "covariates_msd.xlsx"), col_names = TRUE)
msd.exerome.dat <- exerome.dat %>%
  dplyr::left_join(clinical_data, by = "subject") %>%
  dplyr::mutate(t.factor = factor(time))
save(msd.exerome.dat, file = here("data", "msd.exerome.dat.rda"))

# ---- linear mixed model ------------------------------------------------------
registerDoParallel(cores = max(1, parallel::detectCores() - 1))
res.msd.linear <- mixed_anova_parallel(
  msd.exerome.dat, "t.factor", unique(prot.label$Assay), "+ (1|subject)",
  covariates = c("age", "vo2_max")) %>%
  dplyr::left_join(prot.label, by = c("outcome" = "Assay")) %>%
  dplyr::rename(Assay = outcome)
stopImplicitCluster()
res.msd.linear <- as.data.table(res.msd.linear)

save(res.msd.linear, file = here("data", "res.msd.linear.rda"))
cat("Analysis 01 done: wrote metabolic.cytokine.prot.label, msd.exerome.dat, res.msd.linear (",
    nrow(res.msd.linear), "analytes )\n")
