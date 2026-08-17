# =============================================================================
# Analysis 08 — Plasma Olink HPA tissue / cell-type enrichment per cluster.
#
# Joins the plasma LMM/cluster results to Human Protein Atlas (HPA) tissue and
# single-cell annotation, then Fisher-tests each temporal cluster for enrichment
# of tissue- and cell-type-specific proteins.
#
# Input:  data/res.olink.linear.rda (analysis 04), data-raw/hpa_24.tsv
# Output: data/olink_exerome_hpa.rda, data/sig_vs_celltype_olink.rda,
#         data/sig_vs_tissue_cnt_olink.rda, data/sig_vs_celltype_cnt_olink.rda,
#         data/cluster_vs_tissue_cnt_olink.rda (Fig 3g),
#         data/cluster_vs_cell_cnt_olink.rda   (Fig 3h),
#         data/top_tissue_per_protein_olink.rda (Fig 3f)
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(stringr); library(data.table) })

load(here("data", "res.olink.linear.rda"))
exerome_res_raw <- res.olink.linear
columns_to_compare <- c("pval.exposure0", "pval.exposure0.5", "pval.exposure1", "pval.exposure3", "pval.exposure24")

exerome_res <- exerome_res_raw %>% rowwise() %>%
  mutate(pval.exposure = columns_to_compare[which.min(c_across(all_of(columns_to_compare)))],
         beta.exposure = str_replace(pval.exposure, "pval", "beta"),
         b.high = case_when(
           beta.exposure == "beta.exposure0" ~ as.character(beta.exposure0),
           beta.exposure == "beta.exposure0.5" ~ as.character(beta.exposure0.5),
           beta.exposure == "beta.exposure1" ~ as.character(beta.exposure1),
           beta.exposure == "beta.exposure3" ~ as.character(beta.exposure3),
           beta.exposure == "beta.exposure24" ~ as.character(beta.exposure24),
           TRUE ~ NA_character_)) %>% ungroup() %>%
  mutate(significant = as.factor(ifelse(fdr.aov < 0.05, 1, 0))) %>%
  select(OlinkID, Assay, UniProt, fdr.aov, cluster, b.high)
exerome_res <- exerome_res %>% arrange(Assay, fdr.aov) %>% filter(!duplicated(Assay)) %>%
  mutate(cluster = factor(replace_na(cluster, 0), levels = c(0:6)))
exerome_res$significant <- ifelse(exerome_res$fdr.aov < 0.05, 1, 0)
exerome_res$cluster <- replace_na(exerome_res$cluster, 0)

hpa_data <- utils::read.csv(here("data-raw", "hpa_24.tsv"), sep = "\t") %>% unique()
mult_uniprot <- hpa_data %>% filter(str_detect(Uniprot, ", ")) %>% separate_rows(Uniprot, sep = ", ")
hpa_data <- bind_rows(hpa_data %>% filter(!str_detect(Uniprot, ", ")), mult_uniprot)

exerome_hpa <- exerome_res %>% left_join(hpa_data, by = c("Assay" = "Gene"), multiple = "all")
exerome_hpa$RNA.tissue.specific.nTPM <- as.list(strsplit(exerome_hpa$RNA.tissue.specific.nTPM, ";"))
exerome_hpa <- exerome_hpa %>% unnest(cols = RNA.tissue.specific.nTPM) %>%
  separate(RNA.tissue.specific.nTPM, into = c("tissue", "nTPM"), sep = ": ", extra = "merge", fill = "right") %>%
  mutate(tissue = ifelse(is.na(tissue), "non-specific", tissue))
olink_exerome_hpa <- exerome_hpa
save(olink_exerome_hpa, file = here("data", "olink_exerome_hpa.rda"))

# ---- sig_vs_celltype_olink ---------------------------------------------------
sig_vs_celltype <- olink_exerome_hpa %>%
  mutate(RNA.single.cell.type.specific.nTPM = ifelse(is.na(RNA.single.cell.type.specific.nTPM), NA,
                                                     as.character(RNA.single.cell.type.specific.nTPM))) %>%
  mutate(RNA.single.cell.type.specific.nTPM = strsplit(RNA.single.cell.type.specific.nTPM, ";")) %>%
  unnest(cols = RNA.single.cell.type.specific.nTPM) %>%
  separate(RNA.single.cell.type.specific.nTPM, into = c("celltype", "nTPM"), sep = ": ", extra = "merge", fill = "right")
sig_vs_celltype_olink <- sig_vs_celltype
save(sig_vs_celltype_olink, file = here("data", "sig_vs_celltype_olink.rda"))

# ---- Fisher: cluster x tissue -----------------------------------------------
ft <- olink_exerome_hpa %>% mutate(tissue = as.character(tissue)) %>% unnest(tissue) %>% separate_rows(tissue, sep = ";")
cvt <- ft %>% group_by(cluster, tissue) %>% summarise(in_cluster = n(), .groups = "drop")
cvt <- cvt %>% left_join(ft %>% group_by(tissue) %>% summarise(total_in_tissue = n(), .groups = "drop"), by = "tissue") %>%
  mutate(not_in_cluster = total_in_tissue - in_cluster) %>%
  left_join(ft %>% group_by(cluster) %>% summarise(total_in_cluster = n(), .groups = "drop"), by = "cluster") %>%
  mutate(in_cluster_not_tissue = total_in_cluster - in_cluster,
         not_in_cluster_not_tissue = nrow(ft) - (in_cluster + not_in_cluster + in_cluster_not_tissue)) %>%
  rowwise() %>%
  mutate(fisher_p = fisher.test(matrix(c(in_cluster, not_in_cluster, in_cluster_not_tissue, not_in_cluster_not_tissue), ncol = 2, byrow = TRUE), alternative = "greater")$p.value,
         fisher_OR = fisher.test(matrix(c(in_cluster, not_in_cluster, in_cluster_not_tissue, not_in_cluster_not_tissue), ncol = 2, byrow = TRUE), alternative = "greater")$estimate) %>%
  ungroup() %>% group_by(cluster) %>% mutate(padj = p.adjust(fisher_p, method = "fdr")) %>% ungroup()
cluster_vs_tissue_cnt_olink <- cvt
save(cluster_vs_tissue_cnt_olink, file = here("data", "cluster_vs_tissue_cnt_olink.rda"))
sig_vs_tissue_cnt_olink <- cluster_vs_tissue_cnt_olink
save(sig_vs_tissue_cnt_olink, file = here("data", "sig_vs_tissue_cnt_olink.rda"))

# ---- Fisher: cluster x cell type --------------------------------------------
fc <- sig_vs_celltype_olink %>% mutate(celltype = as.character(celltype)) %>% unnest(celltype) %>% separate_rows(celltype, sep = ";")
cvc <- fc %>% group_by(cluster, celltype) %>% summarise(in_cluster = n(), .groups = "drop")
cvc <- cvc %>% left_join(fc %>% group_by(celltype) %>% summarise(total_in_cell = n(), .groups = "drop"), by = "celltype") %>%
  mutate(not_in_cluster = total_in_cell - in_cluster) %>%
  left_join(fc %>% group_by(cluster) %>% summarise(total_in_cluster = n(), .groups = "drop"), by = "cluster") %>%
  mutate(in_cluster_not_cell = total_in_cluster - in_cluster,
         not_in_cluster_not_cell = nrow(fc) - (in_cluster + not_in_cluster + in_cluster_not_cell)) %>%
  rowwise() %>%
  mutate(fisher_p = fisher.test(matrix(c(in_cluster, not_in_cluster, in_cluster_not_cell, not_in_cluster_not_cell), ncol = 2, byrow = TRUE), alternative = "greater")$p.value,
         fisher_OR = fisher.test(matrix(c(in_cluster, not_in_cluster, in_cluster_not_cell, not_in_cluster_not_cell), ncol = 2, byrow = TRUE), alternative = "greater")$estimate) %>%
  ungroup() %>% group_by(cluster) %>% mutate(padj = p.adjust(fisher_p, method = "fdr")) %>% ungroup()
cluster_vs_cell_cnt_olink <- cvc
save(cluster_vs_cell_cnt_olink, file = here("data", "cluster_vs_cell_cnt_olink.rda"))
sig_vs_celltype_cnt_olink <- cluster_vs_cell_cnt_olink
save(sig_vs_celltype_cnt_olink, file = here("data", "sig_vs_celltype_cnt_olink.rda"))

# ---- top tissue per protein: label by first significant timepoint, join HPA --
tissue_order <- c('brain','retina','thyroid gland','parathyroid gland','adrenal gland','pituitary gland','lung','olinkry gland','esophagus','tongue','stomach','intestine','liver','gallbladder','pancreas','kidney','urinary bladder','testis','epididymis','seminal vesicle','prostate','vagina','ovary','fallopian tube','endometrium','cervix','placenta','breast','heart muscle','smooth muscle','skeletal muscle','adipose tissue','skin','lymphoid tissue','bone marrow','non-specific')
sig_vs_tissue <- olink_exerome_hpa %>% distinct(Assay, .keep_all = TRUE) %>%
  mutate(tissue = case_when(tissue == 'choroid plexus' ~ 'brain', tissue == 'stomach 1' ~ 'stomach',
                            tissue == 'skin 1' ~ 'skin', tissue == 'endometrium 1' ~ 'endometrium', TRUE ~ tissue))
sig_vs_tissue$significant <- as.character(sig_vs_tissue$significant)
top_tissue_first <- sig_vs_tissue %>% arrange(Assay, desc(nTPM)) %>% filter(!duplicated(Assay)) %>%
  mutate(tissue = replace_na(tissue, "non-specific"))
top_tissue_first$tissue <- factor(top_tissue_first$tissue, levels = tissue_order)
top_tissue_first <- top_tissue_first %>% group_by(tissue) %>% mutate(tissue_rank = row_number()) %>% ungroup() %>%
  mutate(b.high = as.numeric(as.character(b.high)),
         b_high_dir = case_when(is.na(b.high) ~ NA_character_, b.high > 0 ~ "up", b.high < 0 ~ "down", b.high == 0 ~ "neutral"),
         single_tissue = if_else(RNA.tissue.distribution == "Detected in single", 1, 0),
         tissue_enriched = if_else(RNA.tissue.specificity == "Tissue enriched", 1, 0))

new_df <- res.olink.linear %>% rowwise() %>%
  mutate(pval.exposure = columns_to_compare[which.min(c_across(all_of(columns_to_compare)))],
         beta.exposure = str_replace(pval.exposure, "pval", "beta"),
         b.high = case_when(
           beta.exposure == "beta.exposure0" ~ as.character(beta.exposure0),
           beta.exposure == "beta.exposure0.5" ~ as.character(beta.exposure0.5),
           beta.exposure == "beta.exposure1" ~ as.character(beta.exposure1),
           beta.exposure == "beta.exposure3" ~ as.character(beta.exposure3),
           beta.exposure == "beta.exposure24" ~ as.character(beta.exposure24),
           TRUE ~ NA_character_)) %>% ungroup()
selected_cols <- c("ensembl_ids", "Assay", "UniProt", "fdr.aov", "cluster", "b.high", columns_to_compare)
fr <- new_df %>% select(all_of(selected_cols))
und <- fr %>% filter(str_detect(UniProt, "_"))
if (nrow(und) > 0) fr <- fr %>% filter(!str_detect(UniProt, "_")) %>% bind_rows(und %>% separate_rows(UniProt, sep = "_"))
for (col in columns_to_compare) fr[[paste0(col, "_BH")]] <- p.adjust(fr[[col]], method = "BH")
fr <- as.data.table(fr); fr[, first_sig := NA_character_]
for (i in seq_len(nrow(fr))) {
  below <- FALSE
  for (col in columns_to_compare) { bh_col <- paste0(col, "_BH"); pv <- fr[i, ..bh_col]
    if (!is.na(pv) && pv <= 0.05) { fr[i, first_sig := col]; below <- TRUE; break } }
  if (!below && fr[i, fdr.aov] <= 0.05) {
    sm <- min(fr[i, ..columns_to_compare], na.rm = TRUE)
    fr[i, first_sig := names(fr)[which(fr[i, ..columns_to_compare] == sm)]] }
}
fr[, first_sig := sub("pval.exposure", "", first_sig)]; fr[, first_sig := factor(first_sig)]
fr <- fr[order(fr$Assay, fr$fdr.aov)]; fr <- fr[!duplicated(fr$Assay)]
fr <- fr %>% left_join(top_tissue_first, by = "Assay") %>% mutate(tissue = replace_na(tissue, "non-specific"))
top_tissue_per_protein_olink <- fr
save(top_tissue_per_protein_olink, file = here("data", "top_tissue_per_protein_olink.rda"))

cat("Analysis 08 done: plasma Olink HPA tissue/cell enrichment + top-tissue objects saved.\n")
