# =============================================================================
# Source Data — Supplementary Figure 8 (plasma per-participant protein examples).
#
# Input:  data/res.olink.linear.rda, data/olink.exerome.dat.rda
# Output: source_data/SourceData_SupplementaryFigure8.xlsx
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(writexl) })
dir.create(here("source_data"), showWarnings = FALSE)

load(here("data", "res.olink.linear.rda")); load(here("data", "olink.exerome.dat.rda"))
# rename Olink OID columns -> gene (Assay) names
mapping <- res.olink.linear %>% dplyr::distinct(OlinkID, Assay) %>%
  dplyr::filter(OlinkID %in% colnames(olink.exerome.dat))
npx.data <- olink.exerome.dat
mi <- match(colnames(npx.data), mapping$OlinkID)
colnames(npx.data)[!is.na(mi)] <- mapping$Assay[mi[!is.na(mi)]]
npx.data <- npx.data[, !duplicated(colnames(npx.data))]
npx.data$time    <- as.numeric(as.character(npx.data$time_label))
npx.data$subject <- npx.data$subject

# a: two representative proteins (lowest q) from each cluster 1-6
examples <- res.olink.linear %>% dplyr::filter(!is.na(cluster)) %>%
  dplyr::group_by(cluster) %>% dplyr::arrange(fdr.aov) %>% dplyr::slice_head(n = 2) %>%
  dplyr::ungroup() %>% dplyr::arrange(cluster) %>% dplyr::pull(Assay)

src_long <- function(prots) npx.data %>%
  dplyr::select(time, subject, dplyr::all_of(intersect(prots, colnames(npx.data)))) %>%
  dplyr::filter(time <= 3) %>%
  tidyr::pivot_longer(-c(time, subject), names_to = "protein", values_to = "npx") %>% as.data.frame()

writexl::write_xlsx(list(a_protein_examples = src_long(examples),
                         b_GPCRs = src_long(c("GPR37", "GPR158"))),
                    here("source_data", "SourceData_SupplementaryFigure8.xlsx"))
cat("Source Data Supplementary Figure 8 written: a_protein_examples, b_GPCRs\n")
