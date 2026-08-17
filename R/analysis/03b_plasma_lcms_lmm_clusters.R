# =============================================================================
# Plasma LC-MS/MS — linear mixed model + temporal clusters (Discovery cohort).
#
# Baseline (t = -1) z-normalisation of the plasma LC-MS proteome, then a linear
# mixed model per protein (random subject intercept; age + VO2max covariates).
# Proteins with a time effect at fdr.aov < 0.05 are grouped into temporal
# clusters with navmix (K = 4, seed 42).
#
# Input:  data/plasma_npx_data.rda, data/prot.label.plasma.rda, data-raw/clinical_data.xlsx
# Output: data/exerome.dat.plasma.rda, data/prot.label.plasma.rda,
#         data/res.plasma.linear.rda
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr)
  library(data.table); library(doParallel); library(navmix)
})
source(here("R/package_loading.R")); source(here("R/functions_loading.R"))

load(here("data", "plasma_npx_data.rda"));   npx.data   <- plasma_npx_data
load(here("data", "prot.label.plasma.rda")); prot.label <- prot.label.plasma
clinical_data <- readxl::read_xlsx(here("data-raw", "clinical_data.xlsx"))
clinical_data$subject <- gsub("^EX", "", clinical_data$subject)

exerome.norm <- npx.data %>%
  dplyr::filter(time == -1) %>%
  tidyr::pivot_longer(cols = -c(time, label), names_to = "Assay", values_to = "value") %>%
  dplyr::group_by(Assay) %>%
  dplyr::summarize(mean.value = mean(value, na.rm = TRUE), sd.value = sd(value, na.rm = TRUE), .groups = "drop")

npx.data <- npx.data %>% dplyr::rename(subject = label)
exerome.dat <- npx.data
for (j in seq_len(nrow(exerome.norm))) {
  exerome.dat <- exerome.dat %>%
    dplyr::mutate(!!exerome.norm$Assay[j] := (!!sym(exerome.norm$Assay[j]) - exerome.norm$mean.value[j]) / exerome.norm$sd.value[j])
}

exerome.dat <- exerome.dat %>%
  dplyr::left_join(clinical_data, by = "subject") %>%
  dplyr::mutate(t.factor = factor(time), time_label = as.numeric(as.character(time)))
exerome.dat.plasma <- exerome.dat
save(exerome.dat.plasma, file = here("data", "exerome.dat.plasma.rda"))

prot.label.plasma <- prot.label
save(prot.label.plasma, file = here("data", "prot.label.plasma.rda"))

registerDoParallel(cores = max(1, detectCores() - 1))
res.plasma.linear <- mixed_anova_parallel(exerome.dat, "t.factor", unique(prot.label$Assay),
                                          "+ (1|subject)", covariates = c("age", "vo2_max"))
stopImplicitCluster()

res.plasma.linear <- res.plasma.linear %>%
  left_join(prot.label, by = c("outcome" = "Assay")) %>%
  mutate(fdr.aov = p.adjust(pval.aov.t.factor, method = "BH"),
         fdr.age = p.adjust(pval.aov.age, method = "BH"),
         fdr.vo2_max = p.adjust(pval.aov.vo2_max, method = "BH"))
res.plasma.linear <- as.data.table(res.plasma.linear) %>%
  dplyr::rename(Assay = outcome) %>%
  dplyr::distinct(Assay, .keep_all = TRUE)

cl.data <- res.plasma.linear %>%
  dplyr::filter(fdr.aov < .05) %>%
  dplyr::select(Assay, matches("beta|se|pval")) %>%
  as.data.frame()

if (nrow(cl.data) > 0) {
  rownames(cl.data) <- cl.data$Assay
  tps <- c(0, 0.5, 1, 3, 24)
  for (j in tps) {
    cl.data[, paste0("zscore.", j)] <- cl.data[, paste0("beta.exposure", j)] / cl.data[, paste0("se.exposure", j)]
  }
  set.seed(42)
  fit <- navmix(as.matrix(dplyr::select(cl.data, starts_with("zscore."))),
                K = 4, plot = FALSE, select_K = FALSE, reorder_traits = FALSE)
  cl.data$cluster <- fit$fit$z

  # Optional relabel to match published cluster IDs if sheet A is available.
  tmp.cl <- data.frame(Assay = rownames(cl.data), cluster = as.integer(cl.data$cluster))
  pub_try <- try(suppressWarnings(readxl::read_xlsx(here("data-raw", "supplementary_table_2.xlsx"), sheet = "A")), silent = TRUE)
  if (!inherits(pub_try, "try-error") && all(c("Assay", "cluster") %in% colnames(pub_try))) {
    pub_cl <- setNames(suppressWarnings(as.integer(pub_try$cluster)), pub_try$Assay)
    ov <- table(cl.data$cluster, pub_cl[rownames(cl.data)])
    if (ncol(ov) > 0 && nrow(ov) > 0) {
      relabel <- setNames(as.integer(colnames(ov)[max.col(ov, ties.method = "first")]), rownames(ov))
      tmp.cl <- data.frame(Assay = rownames(cl.data), cluster = unname(relabel[as.character(cl.data$cluster)]))
    }
  }
  res.plasma.linear <- res.plasma.linear %>% left_join(tmp.cl, by = "Assay")
} else {
  res.plasma.linear$cluster <- NA_integer_
}

res.plasma.linear <- res.plasma.linear %>% arrange(fdr.aov)
save(res.plasma.linear, file = here("data", "res.plasma.linear.rda"))

cat(sprintf("Plasma LC-MS: %d proteins, %d exercise-regulated (fdr<0.05). Cluster sizes:\n",
            nrow(res.plasma.linear), sum(res.plasma.linear$fdr.aov < 0.05, na.rm = TRUE)))
print(table(res.plasma.linear$cluster))
