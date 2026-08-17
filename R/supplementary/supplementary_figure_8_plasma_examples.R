# =============================================================================
# Supplementary Figure 8 — Examples of specific proteins with robust responses across participants for each plasma proteome dynamic cluster (related to Fig. 3).
#
#   a  per-participant NPX for two representative proteins from each cluster (1-6)
#   b  two G-protein-coupled receptors regulated by acute exercise (GPR37, GPR158)
#
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2); library(Polychrome); library(gridExtra)
})
OUT <- here("figures", "supplementary"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

load(here("data", "res.olink.linear.rda"))
load(here("data", "olink.exerome.dat.rda"))

# rename Olink OID columns -> gene (Assay) names
mapping <- res.olink.linear %>% distinct(OlinkID, Assay) %>% filter(OlinkID %in% colnames(olink.exerome.dat))
npx.data <- olink.exerome.dat
mi <- match(colnames(npx.data), mapping$OlinkID)
colnames(npx.data)[!is.na(mi)] <- mapping$Assay[mi[!is.na(mi)]]
npx.data <- npx.data[, !duplicated(colnames(npx.data))]
npx.data$time  <- as.numeric(as.character(npx.data$time_label))
npx.data$label <- npx.data$subject                        # participant id

participant_levels <- sort(unique(npx.data$label))
my_colors <- setNames(palette36.colors(length(participant_levels)), participant_levels)

# ---- original per-protein plotting ------------------------------------------
preprocess_data <- function(data, protein) {
  data %>% dplyr::select(time, label, all_of(protein)) %>%
    dplyr::rename(protein_response = !!sym(protein)) %>%
    dplyr::mutate(time = as.numeric(as.character(time))) %>% tidyr::drop_na()
}
# one protein panel, sharing the y-limits passed in (as in the original)
protein_panel <- function(protein, ylim) {
  d <- preprocess_data(npx.data, protein)
  ggplot(d, aes(x = time, y = protein_response, color = label, group = label)) +
    geom_line(size = 0.5, alpha = 0.7) +
    labs(title = protein, x = "Time [h]", y = "npx") +
    theme_bw() +
    scale_color_manual(values = my_colors) +
    ylim(ylim[1], ylim[2]) + xlim(-1, 3) +
    theme(legend.title = element_blank(), panel.grid = element_blank(), legend.position = "none",
          text = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1))
}
# build a panel for each protein, sharing y-limits within its pair
panels_for <- function(proteins) {
  present <- intersect(proteins, colnames(npx.data))
  out <- list()
  for (i in seq(1, length(present), by = 2)) {
    pr <- present[i:min(i + 1, length(present))]
    vals <- unlist(lapply(pr, function(p) preprocess_data(npx.data, p)$protein_response))
    yl <- range(vals, na.rm = TRUE)
    for (p in pr) out[[p]] <- protein_panel(p, yl)
  }
  out
}

# a: two representative proteins (lowest q) from each cluster 1-6
examples <- res.olink.linear %>% filter(!is.na(cluster)) %>%
  group_by(cluster) %>% arrange(fdr.aov) %>% slice_head(n = 2) %>% ungroup() %>% arrange(cluster) %>% pull(Assay)
pa <- panels_for(examples)
ggsave(file.path(OUT, "S8a_plasma_protein_examples.pdf"),
       arrangeGrob(grobs = pa, ncol = 4), width = 16.5, height = 12, units = "cm", dpi = 600)

# b: two exercise-regulated GPCRs
pb <- panels_for(c("GPR37", "GPR158"))
ggsave(file.path(OUT, "S8b_plasma_gpcrs.pdf"),
       arrangeGrob(grobs = pb, ncol = 2), width = 5.5, height = 5, units = "cm", dpi = 600)

src_long <- function(prots) npx.data %>%
  dplyr::select(time, subject, dplyr::all_of(intersect(prots, colnames(npx.data)))) %>%
  dplyr::filter(time <= 3) %>%
  tidyr::pivot_longer(-c(time, subject), names_to = "protein", values_to = "npx")
# source data now built by R/source_data/ (per-figure script)
cat("Supplementary Figure 8 (plasma examples) written to figures/supplementary/S8_*.pdf\n")
