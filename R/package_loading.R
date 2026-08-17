suppressPackageStartupMessages({
	library(OlinkAnalyze)
	library(ComplexHeatmap)
	library(circlize)
	library(ggvenn)
	library(multiWGCNA)
	library(clusterProfiler)
	library(ggpubfigs)
	library(WGCNA)
	library(arrow)
	library(dplyr)
	library(tidyr)
	library(stringr)
	library(limma)
	library(tibble)
	library(ggplot2)
	library(ggrepel)
	library(snakecase)
	library(AnnotationDbi)
	library(org.Hs.eg.db)
	library(biomaRt)
	library(data.table)
	library(lmerTest)
	library(igraph)
	library(doMC)
	library(colorspace)
	library(navmix)
	library(gprofiler2)
	library(purrr)
	library(rlang)
	library(patchwork)
	library(tidyverse)
	library(readxl)
	library(scales)
	library(janitor)
	library(doParallel)
	library(foreach)
	library(gtsummary)
	library(ggsci)
	library(devtools)
	library(MOFA2)
	library(pheatmap)
	library(reshape2)
	library(cowplot)
	library(magrittr)
	library(ggpubr)
	library(PhosR)
	library(SummarizedExperiment)
	library(broom)
	library(broom.mixed)
	library(ggthemes)
	library(knitr)
	library(edgeR)
	library(UpSetR)
	library(splines)
	library(here)
	library(RColorBrewer)
	library(ggpmisc)
	library(emo)
})

#install.packages("devtools")
# # Install Bioconductor packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install(c(
#     "WGCNA",
#     "limma",
#     "AnnotationDbi",
#     "org.Hs.eg.db",
#     "biomaRt",
#     "MOFA2",
#     "PhosR",
#     "SummarizedExperiment",
#     "edgeR"
# ))
#
# # Install CRAN packages
# install.packages(c(
#     "OlinkAnalyze", "arrow", "dplyr", "tidyr", "stringr", "tibble", "ggplot2",
#     "ggrepel", "snakecase", "data.table", "lmerTest", "igraph", "doMC",
#     "colorspace", "gprofiler2", "purrr", "rlang", "patchwork", "tidyverse",
#     "readxl", "scales", "pheatmap", "reshape2", "cowplot", "magrittr",
#     "ggpubr", "broom", "broom.mixed", "ggthemes", "UpSetR", "here",
#     "RColorBrewer", "ggpmisc"
# ))

