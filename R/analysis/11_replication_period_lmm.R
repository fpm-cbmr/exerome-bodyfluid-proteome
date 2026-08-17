# =============================================================================
# Replication cohort — linear mixed models with exercise period as a covariate.
# Fits the two-way (group x time), group-specific and three-way (sex x time x
# mode) models over the replication proteome, adjusting for the randomised
# exercise-mode order (period). Writes the fitted tables to results/period/.
# Input:  data/validation.exerome.dat.rda, data/validation_prot.label.rda,
#         data-raw/validation_exercise_mode/exercise_order.xlsx
# =============================================================================
source(here::here("R/package_loading.R"))
source(here::here("R/functions_loading.R"))
suppressMessages({library(here); library(dplyr); library(tidyr); library(readxl)})

OUT <- here("results/period")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

load(here("data/validation.exerome.dat.rda")); exerome.dat <- validation.exerome.dat
load(here("data/validation_prot.label.rda"))
prot.label <- validation_prot.label
proteins <- unique(prot.label$Assay)

# ---- build period covariate from randomization schedule ---------------------
mod_map <- c(CONT = "A", HIIT = "H", RES = "S")
ord <- read_excel(here("data-raw/validation_exercise_mode/exercise_order.xlsx"))
period_lookup <- ord[, c("ID","EXP1","EXP2","EXP3")] %>%
  mutate(label = sub("^ID","",ID)) %>%
  pivot_longer(c(EXP1,EXP2,EXP3), names_to="session", values_to="mod") %>%
  mutate(group = unname(mod_map[trimws(mod)]),
         period = factor(readr::parse_number(session), levels=c(1,2,3))) %>%
  select(label, group, period)
exerome.dat <- exerome.dat %>% mutate(label=as.character(label)) %>%
  left_join(period_lookup, by=c("label","group"))
stopifnot(!any(is.na(exerome.dat$period)))
exerome.dat$period <- factor(exerome.dat$period)
cat("period merged; distribution:\n"); print(table(exerome.dat$period))

covs <- c("age","sex","total_fat_mass","period")
num_cores <- max(1, parallel::detectCores()-1)
doParallel::registerDoParallel(cores = num_cores)

# ---- (1) two-way group x time + period  (Table 2f) --------------------------
t0 <- Sys.time()
res.validation.linear.period <- suppressMessages(mixed_anova_parallel_group(
  dat=exerome.dat, expo="t.factor", outc=proteins, formel="+ (1|label)",
  group_var="group", covariates=covs)) %>%
  mutate(pval.aov.t.factor.adj = p.adjust(pval.aov.t.factor,"BH"),
         pval.aov.group.adj    = p.adjust(pval.aov.group,"BH"))
if ("pval.aov.period" %in% colnames(res.validation.linear.period))
  res.validation.linear.period$pval.aov.period.adj <- p.adjust(res.validation.linear.period$pval.aov.period,"BH")
saveRDS(res.validation.linear.period, file.path(OUT,"res.validation.linear.period.rds"))
write.csv(res.validation.linear.period, file.path(OUT,"SuppTable2f_replication_twoway_period.csv"), row.names=FALSE)
cat(sprintf("two-way+period done (%.1f min)\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))

# ---- (2) mode-stratified (A/H/S) + period  (Table 2g / Fig 5f,5g) -----------
t0 <- Sys.time()
res_group <- bind_rows(
  run_model_for_group(exerome.dat, "A", proteins, "+ (1|label)", covs) %>% mutate(group="A"),
  run_model_for_group(exerome.dat, "H", proteins, "+ (1|label)", covs) %>% mutate(group="H"),
  run_model_for_group(exerome.dat, "S", proteins, "+ (1|label)", covs) %>% mutate(group="S")) %>%
  group_by(group) %>%
  mutate(pval.exposure0.adj = p.adjust(pval.exposure0,"BH"),
         pval.exposure0.5.adj = p.adjust(pval.exposure0.5,"BH"),
         pval.exposure24.adj = p.adjust(pval.exposure24,"BH")) %>% ungroup()
saveRDS(res_group, file.path(OUT,"res_group_specific.period.rds"))
write.csv(res_group, file.path(OUT,"SuppTable2g_replication_modespecific_period.csv"), row.names=FALSE)
cat(sprintf("mode-stratified+period done (%.1f min)\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))

# ---- (3) three-way sex x group x time + period  (Table 2h / ED Fig 2) -------
t0 <- Sys.time()
res.three_way.period <- suppressMessages(mixed_anova_three_way_parallel_group(
  dat=exerome.dat, expo="t.factor", outc=proteins, formel="+ (1|label)",
  group_var="group", covariates=covs)) %>%
  mutate(pval.aov.t.factor.adj = p.adjust(pval.aov.t.factor,"BH"),
         pval.aov.group.adj    = p.adjust(pval.aov.group,"BH"),
         pval.aov.sex.time.adj = p.adjust(`pval.aov.t.factor:sex`,"BH"))
saveRDS(res.three_way.period, file.path(OUT,"res.validation.three_way.period.rds"))
write.csv(res.three_way.period, file.path(OUT,"SuppTable2h_replication_threeway_period.csv"), row.names=FALSE)
cat(sprintf("three-way+period done (%.1f min)\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))

doParallel::stopImplicitCluster()
cat(sprintf("Replication period LMMs done -> results/period/. Sex x time (FDR<=0.05): %d\n",
            sum(res.three_way.period$pval.aov.sex.time.adj <= 0.05, na.rm = TRUE)))
