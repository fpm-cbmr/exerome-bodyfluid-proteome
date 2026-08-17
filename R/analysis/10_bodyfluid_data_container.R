# =============================================================================
# Analysis 10 — Combined cross-fluid data container (se_bodyfluid).
#
# Assembles one gene_fluid x sample matrix by stacking the saliva, urine and
# plasma matrices (rows = "<gene>_<fluid>", columns = "<subject>_<timepoint>"),
# from data/exerome.dat.saliva, data/exerome.dat.urine and data/olink.exerome.dat.
# Fig 4 c/f use cross-fluid Spearman correlations (rank-based).
#
# Input:  data/exerome.dat.saliva.rda, data/prot.label.saliva.rda,
#         data/exerome.dat.urine.rda,  data/prot.label.urine.rda,
#         data/olink.exerome.dat.rda,  data/olink.prot.label.rda
# Output: data/se_bodyfluid.rda  (Fig 4 c/f)
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tibble); library(SummarizedExperiment)
})

# numeric time_label -> the sample-name timepoint token used across fluids
tl_map <- c("-1" = "Pre", "0" = "T0", "0.5" = "T30", "1" = "T60", "3" = "T180", "24" = "T24")

# exerome.dat (samples x proteins) -> data.frame of gene_fluid x "<subject>_<tp>"
to_fluid_mat <- function(fd, protein_cols, fluid, id2gene = NULL) {
  protein_cols <- intersect(protein_cols, names(fd))
  samp <- paste0(fd$subject, "_", tl_map[as.character(fd$time_label)])
  m <- t(as.matrix(fd[, protein_cols, drop = FALSE]))       # proteins x samples
  colnames(m) <- samp
  genes <- if (is.null(id2gene)) rownames(m) else unname(id2gene[rownames(m)])
  keep <- !is.na(genes) & genes != ""
  m <- m[keep, , drop = FALSE]; genes <- genes[keep]
  m <- m[!duplicated(genes), , drop = FALSE]; genes <- genes[!duplicated(genes)]
  rownames(m) <- paste0(genes, "_", fluid)
  as.data.frame(m)
}

# ---- saliva ------------------------------------------------------------------
load(here("data", "exerome.dat.saliva.rda")); load(here("data", "prot.label.saliva.rda"))
saliva_mat <- to_fluid_mat(exerome.dat.saliva, unique(prot.label.saliva$Assay), "saliva")

# ---- urine -------------------------------------------------------------------
load(here("data", "exerome.dat.urine.rda")); load(here("data", "prot.label.urine.rda"))
urine_mat <- to_fluid_mat(exerome.dat.urine, unique(prot.label.urine$Assay), "urine")

# ---- plasma (OlinkID columns -> gene) ----------------------------------------
load(here("data", "olink.exerome.dat.rda")); load(here("data", "olink.prot.label.rda"))
id2gene <- setNames(as.character(olink.prot.label$Assay), as.character(olink.prot.label$OlinkID))
plasma_mat <- to_fluid_mat(olink.exerome.dat, grep("^OID", names(olink.exerome.dat), value = TRUE),
                           "plasma", id2gene = id2gene)

# ---- combine: unify the sample columns, then stack the fluids' genes ---------
all_samples <- sort(unique(c(colnames(saliva_mat), colnames(urine_mat), colnames(plasma_mat))))
align_cols <- function(m) {
  mm <- matrix(NA_real_, nrow(m), length(all_samples), dimnames = list(rownames(m), all_samples))
  mm[, colnames(m)] <- as.matrix(m)
  as.data.frame(mm)
}
combined_data <- rbind(align_cols(saliva_mat), align_cols(urine_mat), align_cols(plasma_mat))

# ---- SummarizedExperiment ----------------------------------------------------
ordered_metadata <- DataFrame(
  sample  = colnames(combined_data),
  subject = sub("_.*$", "", colnames(combined_data)),
  time    = sub("^[^_]+_", "", colnames(combined_data)),
  row.names = colnames(combined_data))
se_bodyfluid <- SummarizedExperiment(assays = list(exprs = as.matrix(combined_data)),
                                           colData = ordered_metadata)
save(se_bodyfluid, file = here("data", "se_bodyfluid.rda"))

cat(sprintf("Analysis 10 done: se_bodyfluid | %d gene_fluid features x %d samples (saliva %d, urine %d, plasma %d).\n",
            nrow(combined_data), ncol(combined_data),
            nrow(saliva_mat), nrow(urine_mat), nrow(plasma_mat)))
