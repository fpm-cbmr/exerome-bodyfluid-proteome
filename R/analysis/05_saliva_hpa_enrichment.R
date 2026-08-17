# =============================================================================
# Analysis 05 — Saliva HPA tissue / cell-type enrichment per temporal cluster.
#
# Joins the saliva LMM/cluster results to Human Protein Atlas (HPA) tissue and
# single-cell annotation, then Fisher-tests each temporal cluster for enrichment
# of tissue- and cell-type-specific proteins.
#
# Input:  data/res.saliva.linear.rda (analysis 02)
#         data-raw/hpa_24.tsv
# Output: data/saliva_exerome_hpa.rda, data/discovery_saliva_hpa.rda,
#         data/sig_vs_celltype_saliva.rda,
#         data/sig_vs_tissue_cnt_saliva.rda, data/sig_vs_celltype_cnt_saliva.rda,
#         data/cluster_vs_tissue_cnt_saliva.rda  (Fig 2e),
#         data/cluster_vs_cell_cnt_saliva.rda     (Fig 2f),
#         data/top_tissue_per_protein_saliva.rda  (Fig 2c)
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(tidyr); library(stringr); library(data.table) })

load(here("data", "res.saliva.linear.rda"))
exerome_res_raw <- res.saliva.linear

# ---- direction of peak change + tidy columns --------------------------------
columns_to_compare <- c("pval.exposure0", "pval.exposure0.5", "pval.exposure1", "pval.exposure3", "pval.exposure24")
new_df <- exerome_res_raw %>% rowwise() %>%
  mutate(pval.exposure = columns_to_compare[which.min(c_across(all_of(columns_to_compare)))]) %>% ungroup()
new_df <- new_df %>% mutate(beta.exposure = str_replace(pval.exposure, "pval", "beta")) %>%
  mutate(b.high = case_when(
    beta.exposure == "beta.exposure0" ~ as.character(beta.exposure0),
    beta.exposure == "beta.exposure0.5" ~ as.character(beta.exposure0.5),
    beta.exposure == "beta.exposure1" ~ as.character(beta.exposure1),
    beta.exposure == "beta.exposure3" ~ as.character(beta.exposure3),
    beta.exposure == "beta.exposure24" ~ as.character(beta.exposure24),
    TRUE ~ NA_character_))
exerome_res <- new_df %>% select(Assay, uniprotswissprot, protein, fdr.aov, cluster, b.high)

exerome_res$significant <- ifelse(exerome_res$fdr.aov < 0.05, 1, 0)
exerome_res$cluster <- replace_na(exerome_res$cluster, 0)
# take the more significant result when a protein appears under multiple IDs
exerome_res <- exerome_res[order(exerome_res$Assay, exerome_res$fdr.aov), ]
exerome_res <- exerome_res[!duplicated(exerome_res$Assay), ]
exerome_res$cluster <- factor(exerome_res$cluster, levels = c(0, 1, 2, 3, 4, 5, 6, 7))
exerome_res$significant <- as.factor(exerome_res$significant)

# ---- HPA annotation ----------------------------------------------------------
hpa_data <- read.csv(here("data-raw", "hpa_24.tsv"), sep = "\t")
mult_uniprot <- hpa_data %>% dplyr::filter(str_detect(hpa_data$Uniprot, ", "))
hpa_data <- hpa_data %>% dplyr::filter(!str_detect(hpa_data$Uniprot, ", "))
mult_uniprot <- separate_rows(mult_uniprot, Uniprot, sep = ", ")
hpa_data <- rbind(hpa_data, mult_uniprot)
hpa_data <- unique(hpa_data)

saliva_exerome_hpa <- exerome_res %>% left_join(hpa_data, by = c("Assay" = "Gene"), multiple = "all")
save(saliva_exerome_hpa, file = here("data", "saliva_exerome_hpa.rda"))

# ---- discovery_saliva_hpa: unnest RNA tissue-specific nTPM -> tissue ---------
exerome_hpa <- saliva_exerome_hpa
exerome_hpa$RNA.tissue.specific.nTPM <- as.list(strsplit(exerome_hpa$RNA.tissue.specific.nTPM, ";"))
exerome_hpa <- exerome_hpa %>% unnest(cols = RNA.tissue.specific.nTPM) %>%
  separate(RNA.tissue.specific.nTPM, c("tissue", "nTPM"), sep = ": ", extra = "merge", fill = "right")
discovery_saliva_hpa <- exerome_hpa
save(discovery_saliva_hpa, file = here("data", "discovery_saliva_hpa.rda"))

# ---- sig_vs_celltype_saliva: unnest single-cell-type nTPM -> celltype --------
sig_vs_celltype <- saliva_exerome_hpa
sig_vs_celltype <- sig_vs_celltype %>%
  mutate(RNA.single.cell.type.specific.nTPM = ifelse(is.na(RNA.single.cell.type.specific.nTPM), NA,
                                                     as.character(RNA.single.cell.type.specific.nTPM))) %>%
  mutate(RNA.single.cell.type.specific.nTPM = strsplit(RNA.single.cell.type.specific.nTPM, ";"))
sig_vs_celltype <- unnest(sig_vs_celltype, cols = RNA.single.cell.type.specific.nTPM) %>%
  separate(RNA.single.cell.type.specific.nTPM, c("celltype", "nTPM"), sep = ": ")
sig_vs_celltype_saliva <- sig_vs_celltype
save(sig_vs_celltype_saliva, file = here("data", "sig_vs_celltype_saliva.rda"))

# ---- Fisher enrichment: cluster x tissue ------------------------------------
dsh <- discovery_saliva_hpa %>% mutate(tissue = as.character(tissue)) %>%
  unnest(tissue) %>% separate_rows(tissue, sep = ";")
cluster_vs_tissue_cnt <- dsh %>% group_by(cluster, tissue) %>% summarise(in_cluster = n(), .groups = "drop")
total_tissue_counts <- dsh %>% group_by(tissue) %>% summarise(total_in_tissue = n(), .groups = "drop")
cluster_vs_tissue_cnt <- cluster_vs_tissue_cnt %>% left_join(total_tissue_counts, by = "tissue") %>%
  mutate(not_in_cluster = total_in_tissue - in_cluster)
total_cluster_counts <- dsh %>% group_by(cluster) %>% summarise(total_in_cluster = n(), .groups = "drop")
cluster_vs_tissue_cnt <- cluster_vs_tissue_cnt %>% left_join(total_cluster_counts, by = "cluster") %>%
  mutate(in_cluster_not_tissue = total_in_cluster - in_cluster)
total_proteins <- nrow(dsh)
cluster_vs_tissue_cnt <- cluster_vs_tissue_cnt %>%
  mutate(not_in_cluster_not_tissue = total_proteins - (in_cluster + not_in_cluster + in_cluster_not_tissue)) %>%
  rowwise() %>%
  mutate(fisher_p = fisher.test(matrix(c(in_cluster, not_in_cluster, in_cluster_not_tissue, not_in_cluster_not_tissue), ncol = 2, byrow = TRUE), alternative = "greater")$p.value,
         fisher_OR = fisher.test(matrix(c(in_cluster, not_in_cluster, in_cluster_not_tissue, not_in_cluster_not_tissue), ncol = 2, byrow = TRUE), alternative = "greater")$estimate) %>%
  ungroup() %>% group_by(cluster) %>% mutate(padj = p.adjust(fisher_p, method = "fdr")) %>% ungroup()
cluster_vs_tissue_cnt_saliva <- cluster_vs_tissue_cnt
save(cluster_vs_tissue_cnt_saliva, file = here("data", "cluster_vs_tissue_cnt_saliva.rda"))
sig_vs_tissue_cnt_saliva <- cluster_vs_tissue_cnt_saliva
save(sig_vs_tissue_cnt_saliva, file = here("data", "sig_vs_tissue_cnt_saliva.rda"))

# ---- Fisher enrichment: cluster x cell type ---------------------------------
dch <- sig_vs_celltype_saliva %>% mutate(celltype = as.character(celltype)) %>%
  unnest(celltype) %>% separate_rows(celltype, sep = ";")
cluster_vs_cell_cnt <- dch %>% group_by(cluster, celltype) %>% summarise(in_cluster = n(), .groups = "drop")
total_cell_counts <- dch %>% group_by(celltype) %>% summarise(total_in_cell = n(), .groups = "drop")
cluster_vs_cell_cnt <- cluster_vs_cell_cnt %>% left_join(total_cell_counts, by = "celltype") %>%
  mutate(not_in_cluster = total_in_cell - in_cluster)
total_cluster_counts <- dch %>% group_by(cluster) %>% summarise(total_in_cluster = n(), .groups = "drop")
cluster_vs_cell_cnt <- cluster_vs_cell_cnt %>% left_join(total_cluster_counts, by = "cluster") %>%
  mutate(in_cluster_not_cell = total_in_cluster - in_cluster)
total_proteins <- nrow(dch)
cluster_vs_cell_cnt <- cluster_vs_cell_cnt %>%
  mutate(not_in_cluster_not_cell = total_proteins - (in_cluster + not_in_cluster + in_cluster_not_cell)) %>%
  rowwise() %>%
  mutate(fisher_p = fisher.test(matrix(c(in_cluster, not_in_cluster, in_cluster_not_cell, not_in_cluster_not_cell), ncol = 2, byrow = TRUE), alternative = "greater")$p.value,
         fisher_OR = fisher.test(matrix(c(in_cluster, not_in_cluster, in_cluster_not_cell, not_in_cluster_not_cell), ncol = 2, byrow = TRUE), alternative = "greater")$estimate) %>%
  ungroup() %>% group_by(cluster) %>% mutate(padj = p.adjust(fisher_p, method = "fdr")) %>% ungroup()
cluster_vs_cell_cnt_saliva <- cluster_vs_cell_cnt
save(cluster_vs_cell_cnt_saliva, file = here("data", "cluster_vs_cell_cnt_saliva.rda"))
sig_vs_celltype_cnt_saliva <- cluster_vs_cell_cnt_saliva
save(sig_vs_celltype_cnt_saliva, file = here("data", "sig_vs_celltype_cnt_saliva.rda"))

# ---- top tissue per protein: step 1 = one row per protein w/ its top tissue --
sig_vs_tissue <- discovery_saliva_hpa %>% distinct(Assay, .keep_all = TRUE) %>%
  mutate(tissue = case_when(tissue == "choroid plexus" ~ "brain", tissue == "stomach 1" ~ "stomach",
                            tissue == "skin 1" ~ "skin", tissue == "endometrium 1" ~ "endometrium", TRUE ~ tissue))
sig_vs_tissue$significant <- as.character(sig_vs_tissue$significant)
top_tissue_first <- sig_vs_tissue %>% arrange(Assay, desc(nTPM)) %>% filter(!duplicated(Assay)) %>%
  arrange(tissue, fdr.aov) %>% mutate(tissue = replace_na(tissue, "non-specific")) %>%
  group_by(tissue) %>% mutate(tissue_rank = row_number()) %>% ungroup()
top_tissue_first$b.high <- as.numeric(as.character(top_tissue_first$b.high))
top_tissue_first <- top_tissue_first %>%
  mutate(b_high_dir = if_else(sign(b.high) == 1, "up", "down"),
         single_tissue = if_else(RNA.tissue.distribution == "Detected in single", 1, 0),
         tissue_enriched = if_else(RNA.tissue.specificity == "Tissue enriched", 1, 0))

# ---- top tissue per protein: step 2 = join back to ALL proteins (Fig 2c) -----
# label each protein by its first significant timepoint, then join HPA tissue calls
new_df2 <- res.saliva.linear %>% rowwise() %>%
  mutate(pval.exposure = columns_to_compare[which.min(c_across(all_of(columns_to_compare)))],
         beta.exposure = str_replace(pval.exposure, "pval", "beta"),
         b.high = case_when(
           beta.exposure == "beta.exposure0" ~ as.character(beta.exposure0),
           beta.exposure == "beta.exposure0.5" ~ as.character(beta.exposure0.5),
           beta.exposure == "beta.exposure1" ~ as.character(beta.exposure1),
           beta.exposure == "beta.exposure3" ~ as.character(beta.exposure3),
           beta.exposure == "beta.exposure24" ~ as.character(beta.exposure24),
           TRUE ~ NA_character_)) %>% ungroup()
selected_cols <- c("ensembl_ids", "Assay", "uniprotswissprot", "fdr.aov", "cluster", "b.high", columns_to_compare)
exerome_res2 <- new_df2 %>% select(all_of(selected_cols))
und <- exerome_res2 %>% filter(str_detect(uniprotswissprot, "_"))
if (nrow(und) > 0) {
  mult_prots <- und %>% separate_rows(uniprotswissprot, sep = "_")
  exerome_res2 <- exerome_res2 %>% filter(!str_detect(uniprotswissprot, "_")) %>% bind_rows(mult_prots)
}
for (col in columns_to_compare) exerome_res2[[paste0(col, "_BH")]] <- p.adjust(exerome_res2[[col]], method = "BH")
exerome_res2 <- as.data.table(exerome_res2)
exerome_res2[, first_sig := NA_character_]
for (i in seq_len(nrow(exerome_res2))) {
  below <- FALSE
  for (col in columns_to_compare) {
    bh_col <- paste0(col, "_BH"); pval <- exerome_res2[i, ..bh_col]
    if (!is.na(pval) && pval <= 0.05) { exerome_res2[i, first_sig := col]; below <- TRUE; break }
  }
  if (!below && exerome_res2[i, fdr.aov] <= 0.05) {
    sm <- min(exerome_res2[i, ..columns_to_compare], na.rm = TRUE)
    col_name <- names(exerome_res2)[which(exerome_res2[i, ..columns_to_compare] == sm)]
    exerome_res2[i, first_sig := col_name]
  }
}
exerome_res2[, first_sig := sub("pval.exposure", "", first_sig)]
exerome_res2[, first_sig := factor(first_sig)]
exerome_res2 <- exerome_res2[order(exerome_res2$Assay, exerome_res2$fdr.aov), ]
exerome_res2 <- exerome_res2[!duplicated(exerome_res2$Assay), ]
exerome_res2 <- exerome_res2 %>%
  left_join(top_tissue_first, by = "Assay") %>% select(-ends_with(".y")) %>%
  rename_with(~ gsub("\\.x$", "", .), everything())
top_tissue_per_protein_saliva <- exerome_res2 %>% mutate(tissue = ifelse(is.na(tissue), "non-specific", tissue))
save(top_tissue_per_protein_saliva, file = here("data", "top_tissue_per_protein_saliva.rda"))

cat("Analysis 05 done: saliva HPA tissue/cell enrichment + top-tissue objects saved.\n")
