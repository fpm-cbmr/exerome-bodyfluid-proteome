# ==============================================================================
# Three-way (sex x time x mode) model for the Replication cohort.
# y ~ t.factor * group * sex + age + total_fat_mass + (1|label)
# The sex x time effect is the t.factor:sex ANOVA term.
# Fit base and +period; save tables. Writes only to period_sensitivity/output/.
# ==============================================================================
suppressMessages({library(here); library(dplyr); library(tidyr); library(readxl)
                  library(lme4); library(lmerTest); library(parallel)})
OUT <- here("results/period")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

load(here("data/validation.exerome.dat.rda")); dat <- validation.exerome.dat
load(here("data/validation_prot.label.rda")); proteins <- unique(validation_prot.label$Assay)

# period covariate
mod_map <- c(CONT="A",HIIT="H",RES="S")
ord <- read_excel(here("data-raw/validation_exercise_mode/exercise_order.xlsx"))
pl <- ord[,c("ID","EXP1","EXP2","EXP3")] %>% mutate(label=sub("^ID","",ID)) %>%
  pivot_longer(c(EXP1,EXP2,EXP3),names_to="s",values_to="mod") %>%
  mutate(group=unname(mod_map[trimws(mod)]),period=factor(readr::parse_number(s),levels=c(1,2,3))) %>%
  select(label,group,period)
dat <- dat %>% mutate(label=as.character(label)) %>% left_join(pl,by=c("label","group"))
dat$t.factor <- factor(dat$t.factor); dat$group <- factor(dat$group); dat$sex <- factor(dat$sex)

# per-protein corrected fit; return the sex x time (t.factor:sex) ANOVA p-value
fit_one <- function(g, add_period){
  d <- dat[, c(g, "t.factor","group","sex","age","total_fat_mass","period","label")]
  names(d)[1] <- "y"; d <- d[is.finite(d$y), ]
  rhs <- "t.factor * group * sex + age + total_fat_mass"
  if (add_period) rhs <- paste(rhs, "+ period")
  f <- as.formula(paste("y ~", rhs, "+ (1|label)"))
  p <- tryCatch({
    m <- suppressMessages(lmer(f, data=d, REML=TRUE,
           control=lmerControl(check.conv.singular=.makeCC("ignore",tol=1e-4))))
    a <- anova(m)
    c(sex_time = a["t.factor:sex","Pr(>F)"],
      sex      = if("sex" %in% rownames(a)) a["sex","Pr(>F)"] else NA,
      time     = a["t.factor","Pr(>F)"])
  }, error=function(e) c(sex_time=NA,sex=NA,time=NA))
  p
}

run_all <- function(add_period, tag){
  t0 <- Sys.time()
  res <- mclapply(proteins, fit_one, add_period=add_period, mc.cores=max(1,detectCores()-1))
  m <- as.data.frame(do.call(rbind, res)); m$protein <- proteins
  m$sex_time_fdr <- p.adjust(m$sex_time, "BH")
  m$time_na <- is.na(m$time)
  saveRDS(m, file.path(OUT, sprintf("corrected_threeway_%s.rds", tag)))
  write.csv(m, file.path(OUT, sprintf("corrected_threeway_%s.csv", tag)), row.names=FALSE)
  cat(sprintf("%s: fit %.1f min | sex x time FDR<0.05 = %d | time-term NA rows = %d/%d\n",
      tag, as.numeric(difftime(Sys.time(),t0,units="mins")),
      sum(m$sex_time_fdr<0.05,na.rm=TRUE), sum(m$time_na), nrow(m)))
  m
}

cat("Running CORRECTED model (no aliasing)...\n")
base_c <- run_all(FALSE, "base")
per_c  <- run_all(TRUE,  "period")
cat("Corrected sex x time model done -> results/period/ (base + period).\n")
