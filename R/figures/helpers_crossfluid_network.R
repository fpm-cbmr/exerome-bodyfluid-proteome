# =============================================================================
# Shared helpers for the cross-fluid Spearman correlation networks
# (Figure 4 panels d/e and Extended Data 1 panels a-d).
#
# Consumes the ORA-aware graphml + ORA-significance flag files produced by
# 14_crossfluid_correlation_networks.R:
#     data/network/network_crossfluid_padj_ORAaware_<tp>.graphml
#     data/network/ORA_sig_crossfluid_padj_<tp>_Module_<m>.csv
# Node colour/shape/size/label and edge weight/width are baked into the graphml;
# module blob colours are recomputed with module_palette() (seed 123) and the FR
# layout is seeded (set.seed(123)), so the panels are deterministic.
# =============================================================================
suppressMessages({ library(igraph); library(ggplot2); library(ggpubr); library(dplyr) })

CF_FLUID_COLORS <- c(plasma = "#4B0082", saliva = "#FFA07A", urine = "#FFD700")
CF_EDGE_POS <- "#0072B2"; CF_EDGE_NEG <- "#D55E00"

cf_graphml_path <- function(tp, net_dir)
  file.path(net_dir, sprintf("network_crossfluid_padj_ORAaware_%s.graphml", tp))

# farthest-point HCL sampling — reproduces the blob colours from the source script
cf_module_palette <- function(k, seed = 123, L_vals = c(70, 80), C_vals = c(55, 70), hue_step = 4) {
  if (k <= 0) return(character(0))
  set.seed(seed)
  hues <- seq(0, 356, by = hue_step)
  cand <- expand.grid(h = hues, L = L_vals, C = C_vals, KEEP.OUT.ATTRS = FALSE)
  cand$hex <- grDevices::hcl(h = cand$h, c = cand$C, l = cand$L)
  cand <- cand[!is.na(cand$hex), , drop = FALSE]
  if (nrow(cand) < k) {
    L2 <- sort(unique(c(L_vals, mean(range(L_vals)) + c(-10, 10))))
    C2 <- sort(unique(c(C_vals, mean(range(C_vals)) + c(-10, 10))))
    h2 <- seq(0, 355, by = max(2, hue_step %/% 2))
    c2 <- expand.grid(h = h2, L = L2, C = C2, KEEP.OUT.ATTRS = FALSE)
    c2$hex <- grDevices::hcl(h = c2$h, c = c2$C, l = c2$L)
    cand <- unique(rbind(cand, c2[!is.na(c2$hex), , drop = FALSE]))
  }
  idx_pool <- seq_len(nrow(cand))
  pool_vivid <- idx_pool[order(-cand$C, -cand$L)][seq_len(min(20, nrow(cand)))]
  sel_idx <- sample(pool_vivid, 1); selected <- cand[sel_idx, c("h", "C", "L"), drop = FALSE]
  dist_to_sel <- function(M, S) {
    dh <- abs(outer(M[, 1], S[, 1], "-")); dh <- pmin(dh, 360 - dh)
    dC <- abs(outer(M[, 2], S[, 2], "-")); dL <- abs(outer(M[, 3], S[, 3], "-"))
    apply(sqrt((1.2 * dh)^2 + (1.0 * dC)^2 + (0.8 * dL)^2), 1, min)
  }
  remaining <- setdiff(idx_pool, sel_idx)
  while (length(sel_idx) < min(k, nrow(cand))) {
    dmin <- dist_to_sel(as.matrix(cand[remaining, c("h", "C", "L")]), as.matrix(selected))
    pick <- remaining[which.max(dmin)]
    sel_idx <- c(sel_idx, pick); selected <- rbind(selected, cand[pick, c("h", "C", "L")])
    remaining <- setdiff(remaining, pick)
  }
  cand$hex[sel_idx][seq_len(k)]
}

# plot one cross-fluid network (ORA-aware blobs, hub labels); returns the graph
cf_plot_network <- function(tp, out_pdf, net_dir, title = paste0("Cross-Fluid Network at ", tp, " h")) {
  g <- read_graph(cf_graphml_path(tp, net_dir), "graphml")
  kept_mods <- sort(unique(V(g)$Module))
  groups <- lapply(kept_mods, function(m) which(V(g)$Module == m)); names(groups) <- paste0("M", kept_mods)
  ora_sig <- sapply(kept_mods, function(m)
    file.exists(file.path(net_dir, sprintf("ORA_sig_crossfluid_padj_%s_Module_%s.csv", tp, m))))
  k_sig <- sum(ora_sig); pal_sig <- cf_module_palette(k_sig)
  mark_cols <- character(length(kept_mods)); mark_border <- character(length(kept_mods))
  si <- which(ora_sig); ni <- which(!ora_sig)
  if (k_sig > 0) { mark_cols[si] <- grDevices::adjustcolor(pal_sig, alpha.f = 0.40); mark_border[si] <- pal_sig }
  if (length(ni)) { mark_cols[ni] <- grDevices::adjustcolor("#B0B0B0", alpha.f = 0.18); mark_border[ni] <- NA }
  # graphml serialises unlabelled nodes as the literal string "NA"; restore real NA
  vlab <- V(g)$label; vlab[vlab == "NA" | vlab == ""] <- NA
  set.seed(123); lay <- layout_with_fr(g, weights = abs(E(g)$weight))
  pdf(out_pdf, width = 12 / 2.54, height = 9 / 2.54); par(mar = c(2, 0.5, 2, 0.5))
  plot(g, layout = lay, vertex.color = V(g)$color, vertex.frame.color = V(g)$frame.color,
       vertex.shape = V(g)$shape, vertex.size = V(g)$size,
       vertex.label = vlab, vertex.label.color = "black", vertex.label.cex = V(g)$label.cex,
       edge.width = E(g)$width, edge.color = ifelse(E(g)$weight > 0, CF_EDGE_POS, CF_EDGE_NEG),
       edge.curved = FALSE, mark.groups = groups, mark.col = mark_cols, mark.border = mark_border,
       mark.expand = 4, margin = c(0, 0, 0, 0))
  title(main = title, cex.main = 0.8, line = 0.2)
  if (k_sig > 0)
    legend("topright", legend = paste0("M", kept_mods[ora_sig]), pch = 15, pt.cex = 1,
           bty = "n", cex = 0.7, col = mark_cols[ora_sig], title = "Sig. ORA modules")
  legend("bottomleft", legend = c("Positive", "Negative"), lty = 1,
         col = c(CF_EDGE_POS, CF_EDGE_NEG), lwd = 2, bty = "n", cex = 0.7, title = "Edge sign")
  legend("bottomright", legend = names(CF_FLUID_COLORS), col = CF_FLUID_COLORS,
         pch = 16, title = "Body Fluid", pt.cex = 0.7, bty = "n", cex = 0.7)
  mtext(paste("Nodes:", vcount(g), "| Edges:", ecount(g), "| Hubs:", sum(V(g)$shape == "square")),
        side = 1, line = 0, cex = 0.7)
  dev.off()
  g
}

# top-3 hub proteins per fluid (degree = # cross-fluid edges) — network inset
cf_top_hub_inset <- function(g, out_pdf) {
  # degree = # edges per node; stable tiebreak (Protein name) so equal-degree
  # proteins resolve deterministically and match run-to-run.
  hub <- data.frame(Protein = sub("_(urine|saliva|plasma)$", "", as.character(V(g)$name)),
                    Fluid = as.character(V(g)$fluid), Edges = as.integer(igraph::degree(g))) %>%
    dplyr::group_by(Fluid) %>% dplyr::arrange(dplyr::desc(Edges), Protein) %>% dplyr::slice(1:3) %>%
    dplyr::ungroup() %>% dplyr::arrange(Edges, Protein) %>%
    dplyr::mutate(Protein = factor(Protein, levels = Protein))
  p <- ggplot(hub, aes(Edges, Protein, fill = Fluid)) +
    geom_col(width = 0.7) + scale_fill_manual(values = CF_FLUID_COLORS) +
    labs(x = "# edges", y = "Top 3 proteins / Fluid") +
    theme_pubr(base_size = 6) + theme(text = element_text(size = 6), legend.position = "none")
  if (exists("theme_strokes")) p <- p + theme_strokes
  ggsave(out_pdf, p, width = 4.5, height = 4, units = "cm", dpi = 600)
  hub
}
