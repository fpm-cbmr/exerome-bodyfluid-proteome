# =============================================================================
# Source Data — Supplementary Figure 6 (a: canonical exerkine z-scores;
# b: MSD-vs-Olink cross-platform cytokine correlation).
#
# Input:  data-raw/exerkine_list.xlsx, data/res.olink.linear.rda,
#         data/msd.exerome.dat.rda, data/olink.exerome.dat.rda
# Output: source_data/SourceData_SupplementaryFigure6.xlsx
# =============================================================================
suppressMessages({ library(here); library(readxl); library(dplyr); library(tidyr); library(tibble); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)
SHEETS <- list()

# ---- a: exerkine z = beta/se per timepoint from the discovery Olink model -----
alias <- c(IL8 = "CXCL8", VEGF = "VEGFA")
to_olink <- function(g) ifelse(g %in% names(alias), alias[g], g)
user_olink <- unname(to_olink(unique(toupper(trimws(readxl::read_excel(here("data-raw", "exerkine_list.xlsx"))$exerkine)))))
mine_olink <- c("IL6","IL10","IL15","CXCL8","LIF","OSM","FGF21","GDF15","TNFSF11",
                "FABP4","GHRL","LEP","VEGFA","DCN","FST","ANGPTL4","BDNF","MSTN","CCL2","SPP1")
exerkine_set <- union(user_olink, mine_olink)   # canonical exerkines to show (from the curated list + literature)
load(here("data", "res.olink.linear.rda"))
res <- as.data.frame(res.olink.linear); res$.rid <- seq_len(nrow(res))
long <- function(cols, val) res %>% dplyr::select(.rid, Assay, dplyr::all_of(cols)) %>%
  tidyr::pivot_longer(-c(.rid, Assay), names_to = "tp", names_pattern = ".*exposure(.*)", values_to = val)
z_tab <- long(grep("^beta\\.exposure", colnames(res), value = TRUE), "beta") %>%
  dplyr::inner_join(long(grep("^se\\.exposure", colnames(res), value = TRUE), "se") %>% dplyr::select(.rid, tp, se), by = c(".rid", "tp")) %>%
  dplyr::inner_join(long(grep("^pval\\.exposure", colnames(res), value = TRUE), "pval") %>% dplyr::select(.rid, tp, pval), by = c(".rid", "tp")) %>%
  dplyr::group_by(tp) %>% dplyr::mutate(fdr = p.adjust(pval, "BH")) %>% dplyr::ungroup() %>% dplyr::mutate(z = beta / se)
SHEETS[["a_exerkine_zscores"]] <- z_tab %>% dplyr::filter(Assay %in% exerkine_set) %>%
  dplyr::mutate(sig = fdr < 0.05) %>%
  dplyr::select(Assay, tp, beta, se, z, pval, fdr, sig) %>%
  dplyr::arrange(Assay, tp) %>% as.data.frame()

# ---- b: MSD vs Olink paired cytokine values + per-cytokine Spearman ----------
load(here("data", "msd.exerome.dat.rda")); load(here("data", "olink.exerome.dat.rda"))
msd.exerome.dat$subject   <- sub("^EX", "", as.character(msd.exerome.dat$subject))
olink.exerome.dat$subject <- sub("^EX", "", as.character(olink.exerome.dat$subject))
shared <- tibble::tribble(~label, ~msd, ~olink_assay,
  "IL-6","IL6","IL6", "IL-8","IL8","CXCL8", "IL-10","IL10","IL10",
  "TNF-a","TNFa","TNF", "IFN-g","IFNy","IFNG", "IL-2","IL2","IL2")
oid_map <- res.olink.linear %>% dplyr::select(OlinkID, Assay) %>% dplyr::distinct() %>%
  dplyr::filter(Assay %in% shared$olink_assay)
msd_long <- msd.exerome.dat %>% dplyr::select(subject, time, dplyr::all_of(shared$msd)) %>%
  tidyr::pivot_longer(-c(subject, time), names_to = "msd", values_to = "msd_value") %>%
  dplyr::left_join(shared, by = "msd")
olink_sel <- olink.exerome.dat %>% dplyr::select(subject, time = time_label, dplyr::all_of(oid_map$OlinkID)) %>%
  tidyr::pivot_longer(-c(subject, time), names_to = "OlinkID", values_to = "olink_value") %>%
  dplyr::left_join(oid_map, by = "OlinkID") %>%
  dplyr::left_join(shared %>% dplyr::select(label, olink_assay), by = c("Assay" = "olink_assay"))
paired <- dplyr::inner_join(msd_long %>% dplyr::select(subject, time, label, msd_value),
                            olink_sel %>% dplyr::select(subject, time, label, olink_value),
                            by = c("subject", "time", "label")) %>%
  dplyr::filter(is.finite(msd_value), is.finite(olink_value))
stats <- paired %>% dplyr::group_by(label) %>%
  dplyr::summarise(n = dplyr::n(), rho = cor(msd_value, olink_value, method = "spearman", use = "complete.obs"),
                   p = suppressWarnings(cor.test(msd_value, olink_value, method = "spearman")$p.value), .groups = "drop") %>%
  dplyr::arrange(desc(rho))
overall <- with(paired, cor.test(msd_value, olink_value, method = "spearman"))
SHEETS[["b_cytokine_correlations"]] <- dplyr::bind_rows(stats,
  tibble::tibble(label = "ALL (pooled)", n = nrow(paired), rho = unname(overall$estimate), p = overall$p.value)) %>% as.data.frame()
SHEETS[["b_paired_values"]] <- as.data.frame(paired)

writexl::write_xlsx(SHEETS, here("source_data", "SourceData_SupplementaryFigure6.xlsx"))
cat("Source Data Supplementary Figure 6 written:", paste(names(SHEETS), collapse = ", "), "\n")
