# =============================================================================
# Source Data — Supplementary Figure 5 (plasma LC-MS clusters a-d + Olink-vs-LC-MS
# correlation e).
#
# Input:  data/res.plasma.linear.rda, data/exerome.dat.plasma.rda,
#         data/plasma_npx_data.rda, data/olink_npx.data.rda, data/res.olink.linear.rda
# Output: source_data/SourceData_SupplementaryFigure5.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

# a-d: per-protein per-timepoint mean z-score for each LC-MS cluster
cluster_profiles <- function(dat, res, id_col = "Assay") {
  res_cl <- res %>% dplyr::filter(!is.na(cluster)) %>% dplyr::select(dplyr::all_of(id_col), cluster)
  ids <- res_cl %>% dplyr::pull(id_col)
  dplyr::bind_rows(lapply(ids, function(m) {
    dat %>% dplyr::group_by(time_label) %>%
      dplyr::summarise(zscore = mean(.data[[m]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(!!id_col := m)
  })) %>% dplyr::left_join(res_cl, by = id_col)
}
load(here("data", "res.plasma.linear.rda")); load(here("data", "exerome.dat.plasma.rda"))
SHEETS[["a_d_lcms_cluster_profiles"]] <- as.data.frame(cluster_profiles(exerome.dat.plasma, res.plasma.linear, "Assay"))

# e: per-protein Olink-vs-LC-MS Spearman correlation for the 170 shared proteins
load(here("data", "plasma_npx_data.rda")); load(here("data", "olink_npx.data.rda")); load(here("data", "res.olink.linear.rda"))
map <- res.olink.linear %>% dplyr::select(OlinkID, Assay) %>% dplyr::distinct()
oid <- colnames(olink_npx.data); ren <- map$Assay[match(oid, map$OlinkID)]; ren[is.na(ren)] <- oid[is.na(ren)]
colnames(olink_npx.data) <- ren
olink_npx.data <- olink_npx.data[, !duplicated(colnames(olink_npx.data))]
common <- setdiff(intersect(colnames(plasma_npx_data), colnames(olink_npx.data)), c("label", "time", "subject"))
pl <- plasma_npx_data %>% dplyr::select(label, time, dplyr::any_of(common)) %>%
  tidyr::pivot_longer(-c(label, time), names_to = "protein", values_to = "ms") %>% dplyr::mutate(time = as.character(time))
ol <- olink_npx.data %>% dplyr::select(label, time, dplyr::any_of(common)) %>%
  tidyr::pivot_longer(-c(label, time), names_to = "protein", values_to = "olink") %>% dplyr::mutate(time = as.character(time))
m <- dplyr::inner_join(pl, ol, by = c("label", "time", "protein")) %>% dplyr::filter(is.finite(ms), is.finite(olink))
ex <- res.olink.linear %>% dplyr::filter(fdr.aov < 0.05) %>% dplyr::pull(Assay) %>% unique()
SHEETS[["e_olink_lcms_correlation"]] <- m %>% dplyr::group_by(protein) %>% dplyr::filter(dplyr::n() >= 4) %>%
  dplyr::summarise(n = dplyr::n(), r = cor(ms, olink, method = "spearman"),
                   r_p = suppressWarnings(cor.test(ms, olink, method = "spearman")$p.value),
                   ms_abundance = mean(ms, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(r_fdr = p.adjust(r_p, "BH"), exercise_reg = protein %in% ex,
                category = dplyr::case_when(exercise_reg & r >= 0 ~ "Exercise-reg, positive r",
                                            exercise_reg & r < 0 ~ "Exercise-reg, negative r",
                                            TRUE ~ "Not exercise-reg")) %>% as.data.frame()

writexl::write_xlsx(SHEETS, here("source_data", "SourceData_SupplementaryFigure5.xlsx"))
cat("Source Data Supplementary Figure 5 written:", paste(names(SHEETS), collapse = ", "), "\n")
