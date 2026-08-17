# =============================================================================
# Urine — linear mixed model + temporal clusters (Discovery cohort).
#
# Baseline (t = -1) z-normalisation of the urine proteome, then a linear mixed
# model per protein (random subject intercept; age + VO2max covariates). Proteins
# with a time effect at fdr.aov < 0.05 are grouped into temporal clusters with
# navmix (K = 4, seed 42).
# The urine cohort has time points 0 and 24 h.
#
# Input:  data/urine_npx_data.rda, data/prot.label.urine.rda, data-raw/clinical_data.xlsx
# Output: data/exerome.dat.urine.rda, data/prot.label.urine.rda, data/res.urine.linear.rda
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr)
  library(data.table); library(doParallel); library(navmix)
})
source(here("R/package_loading.R")); source(here("R/functions_loading.R"))

load(here("data", "urine_npx_data.rda"));   npx.data   <- urine_npx_data
load(here("data", "prot.label.urine.rda")); prot.label <- prot.label.urine
clinical_data <- readxl::read_xlsx(here("data-raw", "clinical_data.xlsx"))
clinical_data$subject <- gsub("^EX", "", clinical_data$subject)

exerome.norm <- npx.data %>% dplyr::filter(time == -1) %>%
  tidyr::pivot_longer(cols = -c(time, label), names_to = "Assay", values_to = "value") %>%
  dplyr::group_by(Assay) %>%
  dplyr::summarize(mean.value = mean(value, na.rm = TRUE), sd.value = sd(value, na.rm = TRUE), .groups = "drop")
npx.data <- npx.data %>% dplyr::rename(subject = label)
exerome.dat <- npx.data
for (j in seq_len(nrow(exerome.norm)))
  exerome.dat <- exerome.dat %>%
    dplyr::mutate(!!exerome.norm$Assay[j] := (!!sym(exerome.norm$Assay[j]) - exerome.norm$mean.value[j]) / exerome.norm$sd.value[j])

exerome.dat <- exerome.dat %>% dplyr::mutate(t.factor = factor(time)) %>% dplyr::left_join(clinical_data, by = "subject")

# figure-ready baseline-normalized frame (numeric time_label for the cluster plotters)
# and the annotation, saved under the names the figure / cross-fluid scripts load.
exerome.dat.urine <- exerome.dat %>%
  dplyr::mutate(time_label = as.numeric(as.character(time))) %>% dplyr::select(-time)
save(exerome.dat.urine, file = here("data", "exerome.dat.urine.rda"))
prot.label.urine <- prot.label
save(prot.label.urine, file = here("data", "prot.label.urine.rda"))

registerDoParallel(cores = max(1, detectCores() - 1))
res.urine.linear <- mixed_anova_parallel(exerome.dat, "t.factor", unique(prot.label$Assay),
                                         "+ (1|subject)", covariates = c("age", "vo2_max"))
stopImplicitCluster()
res.urine.linear <- res.urine.linear %>%
  left_join(prot.label, by = c("outcome" = "Assay")) %>%
  mutate(fdr.aov = p.adjust(pval.aov.t.factor, method = "BH"),
         fdr.age = p.adjust(pval.aov.age, method = "BH"),
         fdr.vo2_max = p.adjust(pval.aov.vo2_max, method = "BH"))
res.urine.linear <- as.data.table(res.urine.linear) %>% dplyr::rename(Assay = outcome) %>%
  dplyr::distinct(Assay, .keep_all = TRUE)

cl.data <- res.urine.linear %>% dplyr::filter(fdr.aov < .05) %>%
  dplyr::select(Assay, matches("beta|se|pval")) %>% as.data.frame()
rownames(cl.data) <- cl.data$Assay
tps <- c(0, 24)
for (j in tps) cl.data[, paste0("zscore.", j)] <- cl.data[, paste0("beta.exposure", j)] / cl.data[, paste0("se.exposure", j)]
set.seed(42)
fit <- navmix(as.matrix(dplyr::select(cl.data, starts_with("zscore."))),
              K = 4, plot = FALSE, select_K = FALSE, reorder_traits = FALSE)
cl.data$cluster <- fit$fit$z

# navmix labels are arbitrary, so give each recomputed cluster the number of the
# published cluster (Supplementary Table 2, sheet C) that its proteins most overlap.
pub_cl <- as.data.frame(suppressWarnings(readxl::read_xlsx(
  here("data-raw", "supplementary_table_2.xlsx"), sheet = "C")))
pub_cl <- setNames(suppressWarnings(as.integer(pub_cl$cluster)), pub_cl$Assay)
ov <- table(cl.data$cluster, pub_cl[rownames(cl.data)])
relabel <- setNames(as.integer(colnames(ov)[max.col(ov, ties.method = "first")]), rownames(ov))
tmp.cl <- data.frame(Assay = rownames(cl.data), cluster = unname(relabel[as.character(cl.data$cluster)]))

res.urine.linear <- res.urine.linear %>% left_join(tmp.cl, by = "Assay") %>% arrange(fdr.aov)
save(res.urine.linear, file = here("data", "res.urine.linear.rda"))
cat(sprintf("Urine: %d proteins, %d exercise-regulated (fdr<0.05). Cluster sizes:\n",
            nrow(res.urine.linear), sum(res.urine.linear$fdr.aov < 0.05, na.rm = TRUE)))
print(table(res.urine.linear$cluster))
