# =============================================================================
# Analysis 02b — Saliva sample t-SNE embedding (Discovery cohort).
#
# Input:  data/exerome.dat.saliva.rda, data/prot.label.saliva.rda
# Output: data/saliva_tsne.rda   (coords + participant labels; used by figure2.R panel a)
# =============================================================================
suppressMessages({ library(here); library(Rtsne) })

load(here("data", "exerome.dat.saliva.rda"))
load(here("data", "prot.label.saliva.rda"))
exerome.dat   <- get(ls(pattern = "exerome.dat")[1])
prot.label <- get(ls(pattern = "prot.label")[1])

assay_cols <- intersect(unique(prot.label$Assay), colnames(exerome.dat))
mat <- as.matrix(exerome.dat[, assay_cols, drop = FALSE])
mat <- mat[, colSums(is.na(mat)) == 0]        # drop proteins with any NA

SEED <- 6
set.seed(SEED)
tsne_out <- Rtsne(mat, dims = 2, perplexity = 30, verbose = FALSE, max_iter = 500)

saliva_tsne <- data.frame(
  tSNE1 = tsne_out$Y[, 1],
  tSNE2 = tsne_out$Y[, 2],
  participant = as.factor(exerome.dat$subject)
)
attr(saliva_tsne, "seed") <- SEED
attr(saliva_tsne, "input_dim") <- dim(mat)
save(saliva_tsne, file = here("data", "saliva_tsne.rda"))
cat(sprintf("Analysis 02b done: saliva_tsne.rda (seed %d, %d samples x %d proteins)\n",
            SEED, nrow(mat), ncol(mat)))
