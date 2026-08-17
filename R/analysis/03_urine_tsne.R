# =============================================================================
# Analysis 03b — Urine sample t-SNE embedding (Discovery cohort).
# Input:  data/exerome.dat.urine.rda, data/prot.label.urine.rda
# Output: data/urine_tsne.rda   (used by figure2.R panel h)
# =============================================================================
suppressMessages({ library(here); library(Rtsne) })
load(here("data", "exerome.dat.urine.rda"))
load(here("data", "prot.label.urine.rda"))
exerome.dat   <- get(ls(pattern = "exerome.dat")[1])
prot.label <- get(ls(pattern = "prot.label")[1])
assay_cols <- intersect(unique(prot.label$Assay), colnames(exerome.dat))
mat <- as.matrix(exerome.dat[, assay_cols, drop = FALSE]); mat <- mat[, colSums(is.na(mat)) == 0]
SEED <- 6; set.seed(SEED)
tsne_out <- Rtsne(mat, dims = 2, perplexity = 20, verbose = FALSE, max_iter = 500)
urine_tsne <- data.frame(tSNE1 = tsne_out$Y[,1], tSNE2 = tsne_out$Y[,2],
                         participant = as.factor(exerome.dat$subject))
attr(urine_tsne, "seed") <- SEED
save(urine_tsne, file = here("data", "urine_tsne.rda"))
cat(sprintf("Analysis 03b done: urine_tsne.rda (seed %d, %d x %d)\n", SEED, nrow(mat), ncol(mat)))
