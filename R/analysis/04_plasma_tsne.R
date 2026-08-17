# =============================================================================
# Analysis 04b — Plasma Olink sample t-SNE (Discovery cohort). (seeded).
# Input:  data/olink.exerome.dat.rda, data/res.olink.linear.rda (assay list)
# Output: data/plasma_tsne.rda   (used by figure3.R panel a)
# =============================================================================
suppressMessages({ library(here); library(Rtsne) })
load(here("data", "olink.exerome.dat.rda"))
load(here("data", "res.olink.linear.rda"))
assay_cols <- intersect(unique(res.olink.linear$OlinkID), colnames(olink.exerome.dat))
mat <- as.matrix(olink.exerome.dat[, assay_cols, drop = FALSE]); mat <- mat[, colSums(is.na(mat)) == 0]
SEED <- 6; set.seed(SEED)
tsne_out <- Rtsne(mat, dims = 2, perplexity = 30, verbose = FALSE, max_iter = 500)
plasma_tsne <- data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2],
                          participant = as.factor(olink.exerome.dat$replicate))
attr(plasma_tsne, "seed") <- SEED
save(plasma_tsne, file = here("data", "plasma_tsne.rda"))
cat(sprintf("Analysis 04b: plasma_tsne.rda (seed %d, %d x %d)\n", SEED, nrow(mat), ncol(mat)))
