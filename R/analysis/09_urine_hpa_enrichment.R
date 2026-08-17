# =============================================================================
# Analysis 09 — Urine HPA tissue / cell-type enrichment per temporal cluster.
#
# Input:  data/res.urine.linear.rda (analysis 03), data-raw/hpa_24.tsv
# Output: data/urine_exerome_hpa.rda, data/sig_vs_celltype_urine.rda,
#         data/sig_vs_tissue_cnt_urine.rda, data/sig_vs_celltype_cnt_urine.rda,
#         data/cluster_vs_tissue_cnt_urine.rda,
#         data/cluster_vs_cell_cnt_urine.rda   (Fig 2 urine cell-type panel),
#         data/top_tissue_per_protein_urine.rda
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(stringr) })

load(here("data", "res.urine.linear.rda"))
columns_to_compare <- c("pval.exposure0", "pval.exposure24")   # urine cohort: pre + 0 h + 24 h
exerome_res <- res.urine.linear %>% rowwise() %>%
  mutate(pval.exposure = columns_to_compare[which.min(c_across(all_of(columns_to_compare)))],
         beta.exposure = str_replace(pval.exposure, "pval", "beta"),
         b.high = case_when(
           beta.exposure == "beta.exposure0" ~ as.character(beta.exposure0),
           beta.exposure == "beta.exposure24" ~ as.character(beta.exposure24),
           TRUE ~ NA_character_)) %>% ungroup() %>%
  select(Assay, uniprotswissprot, protein, fdr.aov, cluster, b.high)
exerome_res$significant <- ifelse(exerome_res$fdr.aov < 0.05, 1, 0)
exerome_res$cluster <- replace_na(exerome_res$cluster, 0)
exerome_res <- exerome_res[order(exerome_res$Assay, exerome_res$fdr.aov), ]
exerome_res <- exerome_res[!duplicated(exerome_res$Assay), ]
exerome_res$cluster <- factor(exerome_res$cluster, levels = c(0, 1, 2, 3, 4, 5))
exerome_res$significant <- as.factor(exerome_res$significant)

hpa_data <- utils::read.csv(here("data-raw", "hpa_24.tsv"), sep = "\t")
mult_uniprot <- hpa_data %>% dplyr::filter(str_detect(hpa_data$Uniprot, ", "))
hpa_data <- hpa_data %>% dplyr::filter(!str_detect(hpa_data$Uniprot, ", "))
mult_uniprot <- separate_rows(mult_uniprot, Uniprot, sep = ", ")
hpa_data <- unique(rbind(hpa_data, mult_uniprot))

exerome_hpa <- exerome_res %>% left_join(hpa_data, by = c("Assay" = "Gene"), multiple = "all")
exerome_hpa$RNA.tissue.specific.nTPM <- as.list(strsplit(exerome_hpa$RNA.tissue.specific.nTPM, ";"))
exerome_hpa <- exerome_hpa %>% unnest(cols = RNA.tissue.specific.nTPM) %>%
  separate(RNA.tissue.specific.nTPM, c("tissue", "nTPM"), sep = ": ", extra = "merge", fill = "right") %>%
  mutate(tissue = ifelse(is.na(tissue), "non-specific", tissue))
urine_exerome_hpa <- exerome_hpa
save(urine_exerome_hpa, file = here("data", "urine_exerome_hpa.rda"))

sig_vs_celltype <- urine_exerome_hpa %>%
  mutate(RNA.single.cell.type.specific.nTPM = ifelse(is.na(RNA.single.cell.type.specific.nTPM), NA,
                                                     as.character(RNA.single.cell.type.specific.nTPM))) %>%
  mutate(RNA.single.cell.type.specific.nTPM = strsplit(RNA.single.cell.type.specific.nTPM, ";")) %>%
  unnest(cols = RNA.single.cell.type.specific.nTPM) %>%
  separate(RNA.single.cell.type.specific.nTPM, c("celltype", "nTPM"), sep = ": ", extra = "merge", fill = "right")
sig_vs_celltype_urine <- sig_vs_celltype
save(sig_vs_celltype_urine, file = here("data", "sig_vs_celltype_urine.rda"))

fisher_cluster <- function(df, feat) {
  df <- df %>% mutate(!!feat := as.character(.data[[feat]])) %>% separate_rows(!!feat, sep = ";")
  cvt <- df %>% group_by(cluster, .data[[feat]]) %>% summarise(in_cluster = n(), .groups = "drop")
  names(cvt)[2] <- feat
  cvt %>%
    left_join(df %>% group_by(.data[[feat]]) %>% summarise(total_in_f = n(), .groups = "drop") %>% setNames(c(feat, "total_in_f")), by = feat) %>%
    mutate(not_in_cluster = total_in_f - in_cluster) %>%
    left_join(df %>% group_by(cluster) %>% summarise(total_in_cluster = n(), .groups = "drop"), by = "cluster") %>%
    mutate(in_cluster_not_f = total_in_cluster - in_cluster,
           not_in_cluster_not_f = nrow(df) - (in_cluster + not_in_cluster + in_cluster_not_f)) %>%
    rowwise() %>%
    mutate(fisher_p = fisher.test(matrix(c(in_cluster, not_in_cluster, in_cluster_not_f, not_in_cluster_not_f), ncol = 2, byrow = TRUE), alternative = "greater")$p.value,
           fisher_OR = fisher.test(matrix(c(in_cluster, not_in_cluster, in_cluster_not_f, not_in_cluster_not_f), ncol = 2, byrow = TRUE), alternative = "greater")$estimate) %>%
    ungroup() %>% group_by(cluster) %>% mutate(padj = p.adjust(fisher_p, method = "fdr")) %>% ungroup()
}
cvt <- fisher_cluster(urine_exerome_hpa, "tissue")
names(cvt)[names(cvt) == "total_in_f"] <- "total_in_tissue"
names(cvt)[names(cvt) == "in_cluster_not_f"] <- "in_cluster_not_tissue"
names(cvt)[names(cvt) == "not_in_cluster_not_f"] <- "not_in_cluster_not_tissue"
cluster_vs_tissue_cnt_urine <- cvt
save(cluster_vs_tissue_cnt_urine, file = here("data", "cluster_vs_tissue_cnt_urine.rda"))
sig_vs_tissue_cnt_urine <- cluster_vs_tissue_cnt_urine
save(sig_vs_tissue_cnt_urine, file = here("data", "sig_vs_tissue_cnt_urine.rda"))

cvc <- fisher_cluster(sig_vs_celltype_urine, "celltype")
names(cvc)[names(cvc) == "total_in_f"] <- "total_in_cell"
names(cvc)[names(cvc) == "in_cluster_not_f"] <- "in_cluster_not_cell"
names(cvc)[names(cvc) == "not_in_cluster_not_f"] <- "not_in_cluster_not_cell"
cluster_vs_cell_cnt_urine <- cvc
save(cluster_vs_cell_cnt_urine, file = here("data", "cluster_vs_cell_cnt_urine.rda"))
sig_vs_celltype_cnt_urine <- cluster_vs_cell_cnt_urine
save(sig_vs_celltype_cnt_urine, file = here("data", "sig_vs_celltype_cnt_urine.rda"))

sig_vs_tissue <- urine_exerome_hpa %>%
  distinct(Assay, .keep_all = TRUE) %>%
  mutate(tissue = case_when(
    tissue == "choroid plexus" ~ "brain",
    tissue == "stomach 1" ~ "stomach",
    tissue == "skin 1" ~ "skin",
    tissue == "endometrium 1" ~ "endometrium",
    TRUE ~ tissue
  ))
top_tissue_first <- sig_vs_tissue %>%
  arrange(Assay, desc(as.numeric(nTPM))) %>%
  filter(!duplicated(Assay)) %>%
  mutate(tissue = replace_na(tissue, "non-specific"))

first_sig_df <- res.urine.linear %>%
  rowwise() %>%
  mutate(pval.exposure = columns_to_compare[which.min(c_across(all_of(columns_to_compare)))]) %>%
  ungroup()
for (col in columns_to_compare) {
  first_sig_df[[paste0(col, "_BH")]] <- p.adjust(first_sig_df[[col]], method = "BH")
}
first_sig_df <- first_sig_df %>%
  mutate(first_sig = dplyr::case_when(
    pval.exposure0_BH <= 0.05 ~ "0",
    pval.exposure24_BH <= 0.05 ~ "24",
    fdr.aov <= 0.05 & pval.exposure0 <= pval.exposure24 ~ "0",
    fdr.aov <= 0.05 & pval.exposure24 < pval.exposure0 ~ "24",
    TRUE ~ NA_character_
  )) %>%
  arrange(Assay, fdr.aov) %>%
  distinct(Assay, .keep_all = TRUE) %>%
  select(Assay, first_sig)

top_tissue_per_protein_urine <- top_tissue_first %>%
  left_join(first_sig_df, by = "Assay")
save(top_tissue_per_protein_urine, file = here("data", "top_tissue_per_protein_urine.rda"))

cat("Analysis 09 done: urine HPA tissue/cell enrichment objects saved.\n")
