# =============================================================================
# Replication cohort — assemble the fitted model tables (results/period/, from
# analyses 11 and 12) into the objects the figures load, adding the BH-adjusted
# ANOVA columns. Writes data/res.validation.linear.pilot.rda,
# data/combined.pilot.res.linear.rda, data/res.validation.three_way.rda and
# data/corrected_threeway_period.rds (Fig 5, ED 2/3).
# ==============================================================================
suppressMessages({library(here); library(dplyr)})
DD  <- here("data")
OUT <- here("results", "period")
PS  <- here("results", "period")

# ---- res.validation.linear.pilot (two-way + period) -------------------------
res.validation.linear.pilot <- readRDS(file.path(OUT,"res.validation.linear.period.rds")) %>%   # analysis 11
  mutate(pval.exposure0.adj      = p.adjust(pval.exposure0,"BH"),
         pval.exposure0.5.adj    = p.adjust(pval.exposure0.5,"BH"),
         pval.exposure24.adj     = p.adjust(pval.exposure24,"BH"),
         pval.aov.t.factor.adj   = p.adjust(pval.aov.t.factor,"BH"),
         pval.aov.group.adj      = p.adjust(pval.aov.group,"BH"),
         pval.aov.t.factor.group.adj = p.adjust(`pval.aov.t.factor:group`,"BH"),
         pval.aov.total_fat_mass.adj = p.adjust(pval.aov.total_fat_mass,"BH"),
         pval.aov.sex.adj        = p.adjust(pval.aov.sex,"BH"),
         pval.aov.age.adj        = p.adjust(pval.aov.age,"BH"))
save(res.validation.linear.pilot, file = file.path(DD,"res.validation.linear.pilot.rda"))

# ---- res_combined_group_specific_time_linear_pilot + combined.pilot.res.linear
# Merge the full-model adjusted aov p-values from the two-way model into the group-specific combined object (needed by the 5f heatmap).
merge_cols <- res.validation.linear.pilot %>%
  select(outcome, pval.aov.group.adj, pval.aov.t.factor.adj, pval.aov.t.factor.group.adj) %>%
  distinct(outcome, .keep_all = TRUE)
res_combined_group_specific_time_linear_pilot <- readRDS(file.path(PS,"res_group_specific.period.rds")) %>%
  mutate(pval.aov.time.adj = p.adjust(pval.aov.time,"BH"),
         pval.aov.sex.adj = p.adjust(pval.aov.sex,"BH"),
         pval.aov.age.adj = p.adjust(pval.aov.age,"BH"),
         pval.aov.total_fat_mass.adj = p.adjust(pval.aov.total_fat_mass,"BH")) %>%
  left_join(merge_cols, by = "outcome")
save(res_combined_group_specific_time_linear_pilot,
     file = file.path(DD,"res_combined_group_specific_time_linear_pilot.rda"))
combined.pilot.res.linear <- res_combined_group_specific_time_linear_pilot
save(combined.pilot.res.linear, file = file.path(DD,"combined.pilot.res.linear.rda"))

# ---- res.validation.three_way ------
res.validation.three_way <- readRDS(file.path(PS,"res.validation.three_way.period.rds"))
# ensure the sex:time adjusted column your plots filter on exists
if (!"pval.aov.sex.time.adj" %in% colnames(res.validation.three_way))
  res.validation.three_way$pval.aov.sex.time.adj <- p.adjust(res.validation.three_way[["pval.aov.t.factor:sex"]],"BH")
save(res.validation.three_way, file = file.path(DD,"res.validation.three_way.rda"))

# ---- validation.exerome.dat  -------------
suppressMessages({library(tidyr); library(readxl)})
load(here("data/validation.exerome.dat.rda"))
mod_map <- c(CONT="A",HIIT="H",RES="S")
ord <- read_excel(here("data-raw/validation_exercise_mode/exercise_order.xlsx"))
pl <- ord[,c("ID","EXP1","EXP2","EXP3")] %>% mutate(label=sub("^ID","",ID)) %>%
  pivot_longer(c(EXP1,EXP2,EXP3),names_to="s",values_to="mod") %>%
  mutate(group=unname(mod_map[trimws(mod)]),period=factor(readr::parse_number(s),levels=c(1,2,3))) %>%
  select(label,group,period)
validation.exerome.dat.period <- validation.exerome.dat %>% mutate(label=as.character(label)) %>%
  left_join(pl,by=c("label","group"))
save(validation.exerome.dat.period, file = file.path(DD,"validation.exerome.dat.period.rda"))

cat("Replication objects packaged into data/:\n")
print(list.files(DD))
cat("\nSanity — sex x time proteins (FDR<=0.05) in packaged three-way:",
    sum(res.validation.three_way$pval.aov.sex.time.adj<=0.05, na.rm=TRUE), "\n")
cat("Mode main effect (group.adj<0.05 & t.factor.adj>0.05):",
    sum(res.validation.linear.pilot$pval.aov.group.adj<0.05 &
        res.validation.linear.pilot$pval.aov.t.factor.adj>0.05, na.rm=TRUE), "\n")

# ---- sex x time model into data/ for Extended Data 2 -------------------------
if (file.exists(file.path(PS, "corrected_threeway_period.rds")))
  file.copy(file.path(PS, "corrected_threeway_period.rds"),
            file.path(DD, "corrected_threeway_period.rds"), overwrite = TRUE)
