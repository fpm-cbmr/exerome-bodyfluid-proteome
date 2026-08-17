# =============================================================================
# Shared helpers for the secreted-exerkine cis-pQTL <-> trait/disease networks
# (main Fig. 6 and Extended Data Figs. 3-4).
#
# cn_build_trait_network_full() builds the trait/disease network from a merged
# colocalization table: Louvain communities + Fruchterman-Reingold layout (both
# seeded, set.seed(123)) and deterministic per-module GO:BP enrichGO. Callers
# supply the merged coloc table, the Olink universe, the trait palette and the
# colouring column (parent_term or category).
# =============================================================================
suppressMessages({
  library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(ggplot2); library(ggpubr); library(igraph); library(scales)
  library(clusterProfiler); library(org.Hs.eg.db); library(fields)
})

# p from Open Targets mantissa/exponent parts
cn_p_from_parts <- function(m, e) {
  m <- suppressWarnings(as.numeric(m)); e <- suppressWarnings(as.numeric(e))
  p <- ifelse(is.finite(m) & is.finite(e), m * (10^e), NA_real_)
  pmax(pmin(ifelse(is.finite(p), p, NA_real_), 1), .Machine$double.xmin)
}

CN_CLUSTER_COLORS <- c("0" = "#B0B0B0", "1" = "#1F78B4", "2" = "#33A02C", "3" = "#E31A1C",
                       "4" = "#FDBF6F", "5" = "#FF7F00", "6" = "#6A3D9A", "7" = "#CAB2D6")

# merge coloc results + trait parents/categories + olink cluster/beta + HPA secretome
cn_prepare_coloc_traits <- function(trait_df, coloc_res, res_olink, hpa_clean) {
  trait_map <- trait_df %>%
    dplyr::rename(reported_traits = dplyr::any_of(c("reported_traits", "reported_trait"))) %>%
    dplyr::select(reported_traits, parent_term, category) %>%
    dplyr::mutate(reported_traits = str_trim(as.character(reported_traits))) %>% dplyr::distinct()
  coloc_df <- coloc_res %>%
    dplyr::mutate(left_trait = str_trim(as.character(left_trait))) %>%
    dplyr::inner_join(trait_map, by = c("left_trait" = "reported_traits")) %>%
    dplyr::mutate(right_geneSymbol = str_trim(as.character(right_geneSymbol)))
  olink_cluster <- res_olink %>% dplyr::filter(fdr.aov < 0.05) %>%
    dplyr::mutate(Assay = str_trim(as.character(Assay))) %>%
    dplyr::select(Assay, cluster, beta.exposure0) %>% dplyr::distinct(Assay, .keep_all = TRUE)
  coloc_df %>%
    dplyr::left_join(olink_cluster, by = c("right_geneSymbol" = "Assay")) %>%
    dplyr::left_join(hpa_clean %>% dplyr::select(gene, secretome_function, `RNA tissue specific nTPM`) %>%
                       dplyr::distinct(), by = c("right_geneSymbol" = "gene"))
}

# build + plot one protein<->trait network; returns a list of the graph + panel data
cn_build_trait_network_full <- function(df, olink_df, out_prefix, title_text, trait_palette,
    color_col = "parent_term", use_secreted = TRUE, min_mod_plot = 10L, alpha_adj = 0.05,
    run_fisher = TRUE, fisher_col = "parent_term", top_n_proteins = 10, top_n_traits = 5,
    cluster_color_palette = CN_CLUSTER_COLORS, emit_cluster_bar = TRUE,
    beta_lab = expression(beta ~ "(0h post-exercise)")) {
  if (use_secreted && "secretome_function" %in% names(df))
    df <- df %>% dplyr::filter(!is.na(secretome_function))

  edges_trait <- df %>% dplyr::mutate(
    edge_dir = dplyr::case_when(
      is.finite(left_beta) & is.finite(right_beta) & left_beta * right_beta >= 0 ~ "risk_like",
      is.finite(left_beta) & is.finite(right_beta) & left_beta * right_beta < 0 ~ "protective_like",
      TRUE ~ NA_character_),
    Protein = str_trim(as.character(right_geneSymbol)), Trait = str_trim(as.character(left_trait)),
    TraitGroup = .data[[color_col]], TraitColor = trait_palette[as.character(.data[[color_col]])]) %>%
    dplyr::filter(!is.na(Protein), Protein != "", !is.na(Trait), Trait != "", !is.na(edge_dir))
  edges_group <- edges_trait %>%
    dplyr::group_by(Protein, Trait, TraitGroup, TraitColor, edge_dir) %>%
    dplyr::summarise(n_hits = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(edge_color = ifelse(edge_dir == "risk_like",
             adjustcolor("#D55E00", alpha.f = 0.5), adjustcolor("#0072B2", alpha.f = 0.5)),
           edge_width = rescale(n_hits, to = c(0.5, 5)))

  protein_nodes <- df %>%
    dplyr::transmute(name = right_geneSymbol, type = "protein",
              beta = as.numeric(beta.exposure0), cluster = as.numeric(cluster)) %>%
    dplyr::distinct() %>% dplyr::filter(cluster != 0)
  beta_trim <- quantile(protein_nodes$beta, c(0.01, 0.99), na.rm = TRUE); beta_max <- max(abs(beta_trim))
  beta_pal <- col_numeric(c("#2166AC", "white", "#B2182B"), domain = c(-beta_max, beta_max))
  protein_nodes <- protein_nodes %>% dplyr::mutate(node_color = beta_pal(pmin(pmax(beta, -beta_max), beta_max)))
  trait_nodes <- edges_group %>% dplyr::transmute(name = Trait, type = "trait",
    node_color = dplyr::if_else(is.na(TraitColor) | TraitColor == "", "#BBBBBB", TraitColor)) %>% dplyr::distinct()
  nodes <- dplyr::bind_rows(protein_nodes, trait_nodes) %>% dplyr::distinct(name, .keep_all = TRUE)
  edges_group <- edges_group %>% dplyr::filter(Protein %in% nodes$name & Trait %in% nodes$name)

  g <- igraph::graph_from_data_frame(edges_group %>% dplyr::transmute(from = Protein, to = Trait,
        color = edge_color, width = edge_width), vertices = nodes, directed = FALSE)
  E(g)$color <- edges_group$edge_color; E(g)$width <- edges_group$edge_width
  V(g)$shape <- ifelse(V(g)$type == "protein", "circle", "square")
  V(g)$frame.color <- "#333333"; V(g)$frame.width <- 0.5
  V(g)$color <- nodes$node_color[match(V(g)$name, nodes$name)]
  comps <- igraph::components(g); g <- igraph::induced_subgraph(g, which(comps$membership == which.max(comps$csize)))
  message("  main component: ", vcount(g), " nodes, ", ecount(g), " edges.")

  g_u <- igraph::as.undirected(g, mode = "collapse",
        edge.attr.comb = list(color = "first", width = "sum", .default = "first"))
  lou <- igraph::cluster_louvain(g_u, weights = E(g_u)$width)
  memb <- as.integer(igraph::membership(lou)); names(memb) <- V(g_u)$name
  V(g)$Module <- memb[match(V(g)$name, names(memb))]

  universe <- olink_df %>% dplyr::mutate(Assay = str_trim(as.character(Assay))) %>%
    dplyr::pull(Assay) %>% unique() %>% { .[nzchar(.)] }
  mods <- sort(unique(V(g)$Module)); ora_sig <- logical(length(mods)); names(ora_sig) <- mods
  for (mid in mods) {
    vids <- which(V(g)$Module == mid & V(g)$type == "protein"); if (length(vids) < 3) next
    ora <- tryCatch(enrichGO(gene = unique(V(g)$name[vids]), universe = universe,
      OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "fdr",
      pvalueCutoff = 0.05, qvalueCutoff = 0.2), error = function(e) NULL)
    if (!is.null(ora) && nrow(as.data.frame(ora)) > 0) {
      dfo <- as.data.frame(ora); pv <- suppressWarnings(as.numeric(dfo$p.adjust))
      ora_sig[as.character(mid)] <- any(is.finite(pv) & pv < alpha_adj)
      if (any(is.finite(pv) & pv < alpha_adj))
        write.csv(dfo[is.finite(pv) & pv < alpha_adj, ],
                  paste0(out_prefix, "_ORA_Module_", mid, ".csv"), row.names = FALSE)
    }
  }

  # module shading (computed BEFORE the Fisher plot)
  kept_mods <- as.integer(names(table(V(g)$Module)[table(V(g)$Module) >= min_mod_plot]))
  sig_mods  <- kept_mods[ora_sig[match(kept_mods, as.integer(names(ora_sig)))] %in% TRUE]
  nonsig    <- setdiff(kept_mods, sig_mods)
  module_palette <- function(k, seed = 123) { set.seed(seed); if (k <= 0) return(character(0))
    h <- seq(0, 330, length.out = max(3, k)); hcl(h = h[seq_len(k)], c = 60, l = 75) }
  pal_sig <- module_palette(length(sig_mods))
  mark_cols <- character(length(kept_mods)); mark_border <- character(length(kept_mods))
  if (length(sig_mods)) { mark_cols[match(sig_mods, kept_mods)] <- adjustcolor(pal_sig, alpha.f = 0.4)
    mark_border[match(sig_mods, kept_mods)] <- pal_sig }
  if (length(nonsig)) { mark_cols[match(nonsig, kept_mods)] <- adjustcolor("#B0B0B0", alpha.f = 0.2)
    mark_border[match(nonsig, kept_mods)] <- NA }
  groups_to_mark <- lapply(kept_mods, function(m) which(V(g)$Module == m)); names(groups_to_mark) <- paste0("M", kept_mods)

  fisher_plot <- NULL; fisher_df <- NULL
  if (run_fisher) {
    node_df <- data.frame(name = V(g)$name, Module = V(g)$Module, type = V(g)$type, stringsAsFactors = FALSE)
    tmap <- df %>% dplyr::select(left_trait, !!rlang::sym(fisher_col)) %>% dplyr::distinct()
    node_df$trait_group <- tmap[[fisher_col]][match(node_df$name, tmap$left_trait)]
    node_df$trait_group[is.na(node_df$trait_group) | node_df$type != "trait"] <- "None"
    fo <- list()
    for (mid in unique(node_df$Module)) {
      mn <- node_df[node_df$Module == mid, ]; on <- node_df[node_df$Module != mid, ]
      for (cat in setdiff(unique(node_df$trait_group), "None")) {
        a <- sum(mn$trait_group == cat); b <- sum(mn$trait_group != cat)
        cc <- sum(on$trait_group == cat); d <- sum(on$trait_group != cat)
        if ((a + cc) > 0) { ft <- fisher.test(matrix(c(a, b, cc, d), nrow = 2))
          fo[[length(fo) + 1]] <- data.frame(Module = mid, TraitGroup = cat, a, b, c = cc, d,
            OR = unname(ft$estimate), pval = ft$p.value) } }
    }
    if (length(fo)) {
      fisher_df <- do.call(rbind, fo); fisher_df$p_adj <- p.adjust(fisher_df$pval, method = "fdr")
      write.csv(fisher_df, paste0(out_prefix, "_Fisher_Module_TraitEnrichment.csv"), row.names = FALSE)
      fpd <- fisher_df %>% dplyr::mutate(log10p = -log10(p_adj), log2OR = log2(OR),
        Sig = p_adj < 0.05) %>% dplyr::filter(Sig & OR > 1) %>% dplyr::mutate(Module = factor(Module))
      fisher_plot <- ggplot(fpd, aes(x = Module, y = log2OR, size = log10p, fill = TraitGroup)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
        geom_point(shape = 21, color = "black", alpha = 0.85, stroke = 0.2) +
        scale_fill_manual(values = trait_palette, na.value = "grey70") +
        scale_size_continuous(range = c(1, 5), name = expression(-log[10]("FDR p-value"))) +
        labs(title = "Trait Enrichment per Network Module", x = "Module", y = expression(log[2]("Odds Ratio"))) +
        guides(fill = "none") + theme_pubr(base_size = 6) +
        scale_x_discrete(labels = function(x) paste0("M", x)) +
        theme(text = element_text(size = 6), panel.grid = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
      if (exists("theme_strokes")) fisher_plot <- fisher_plot + theme_strokes
      ggsave(paste0(out_prefix, "_FisherEnrichment_DotPlot.pdf"), fisher_plot, width = 6.5, height = 8, units = "cm")
    }
  }

  deg_vec <- igraph::degree(g)
  top_proteins <- tibble(name = names(deg_vec[V(g)$type == "protein"]),
    degree = as.numeric(deg_vec[V(g)$type == "protein"])) %>% dplyr::arrange(dplyr::desc(degree)) %>% dplyr::slice_head(n = top_n_proteins)
  top_traits <- tibble(name = names(deg_vec[V(g)$type == "trait"]),
    degree = as.numeric(deg_vec[V(g)$type == "trait"])) %>% dplyr::arrange(dplyr::desc(degree)) %>% dplyr::slice_head(n = top_n_traits)
  prot_trait_comp <- edges_group %>% dplyr::filter(Protein %in% top_proteins$name) %>%
    dplyr::group_by(Protein, TraitGroup) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(Protein) %>% dplyr::mutate(total_n = sum(n)) %>% dplyr::ungroup()
  p1 <- ggplot(prot_trait_comp, aes(reorder(Protein, total_n), n, fill = TraitGroup)) +
    geom_col(position = "stack") + coord_flip() +
    scale_fill_manual(values = trait_palette, na.value = "grey70") + theme_pubr(base_size = 6) +
    labs(title = "Top Proteins: composition of connected traits", x = "Protein", y = "Connections") +
    theme(text = element_text(size = 6), legend.position = "none", panel.grid = element_blank())
  if (exists("theme_strokes")) p1 <- p1 + theme_strokes
  ggsave(paste0(out_prefix, "_TopProteins_StackedTraits.pdf"), p1, width = 7, height = 5, units = "cm")
  edge_cl <- edges_group %>% dplyr::left_join(protein_nodes %>% dplyr::select(name, cluster), by = c("Protein" = "name")) %>%
    dplyr::filter(!is.na(cluster), cluster != 0)
  trait_cluster_comp <- edge_cl %>% dplyr::filter(Trait %in% top_traits$name) %>% dplyr::mutate(cluster = as.character(cluster)) %>%
    dplyr::group_by(Trait, cluster) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(Trait) %>% dplyr::mutate(total_n = sum(n)) %>% dplyr::ungroup()
  # the exercise-temporal-cluster bar (panel c of Fig 6 / ED3a-c). Callers whose
  # "cluster" is not the temporal cluster (e.g. the exercise-mode network) pass
  # emit_cluster_bar = FALSE and build their own composition bar instead.
  if (emit_cluster_bar) {
    p2 <- ggplot(trait_cluster_comp, aes(reorder(Trait, total_n), n, fill = cluster)) +
      geom_col(position = "stack") + coord_flip() +
      scale_fill_manual(values = cluster_color_palette, na.value = "grey70") + theme_pubr(base_size = 6) +
      labs(title = "Top Traits: composition of connected protein clusters", x = "Trait",
           y = "# of Exerkine pQTL colocalizations", fill = "Exercise Cluster") +
      theme(text = element_text(size = 6), legend.position = "right", legend.key.size = unit(2, "mm"), panel.grid = element_blank())
    if (exists("theme_strokes")) p2 <- p2 + theme_strokes
    ggsave(paste0(out_prefix, "_TopTraits_StackedClusters.pdf"), p2, width = 7, height = 5, units = "cm")
  }

  set.seed(123); lay <- igraph::layout_with_fr(g, niter = 1000, grid = "nogrid") * 4
  pdf(paste0(out_prefix, "_network.pdf"), width = 15 / 2.54, height = 15 / 2.54); par(mar = c(2, 0.5, 2, 0.5))
  plot(g, layout = lay, vertex.color = V(g)$color,
       vertex.shape = ifelse(V(g)$type == "protein", "circle", "square"),
       vertex.size = ifelse(V(g)$type == "protein", 1.25, 1.75), vertex.label = NA,
       vertex.frame.color = V(g)$frame.color, edge.width = E(g)$width, edge.color = E(g)$color,
       edge.arrow.size = 0, mark.groups = groups_to_mark, mark.col = mark_cols,
       mark.border = mark_border, mark.expand = 5)
  title(main = title_text, cex.main = 0.9, line = 0.2)
  fields::image.plot(legend.only = TRUE, zlim = c(-beta_max, beta_max),
    col = beta_pal(seq(-beta_max, beta_max, length.out = 100)), legend.shrink = 0.2,
    legend.width = 0.3, horizontal = FALSE, legend.lab = beta_lab)
  legend("bottomright", legend = c("Risk-like", "Protective-like"),
    col = c(adjustcolor("#D55E00", alpha.f = 0.5), adjustcolor("#0072B2", alpha.f = 0.5)),
    lwd = 2, bty = "n", cex = 0.7, title = "Association")
  used_traits <- intersect(unique(na.omit(df[[color_col]])), names(trait_palette))
  if (length(used_traits)) legend("bottomleft", legend = used_traits, col = trait_palette[used_traits],
    pch = 15, bty = "n", cex = 0.6, ncol = 2, title = "Trait categories")
  legend("bottom", legend = c("Protein", "Trait"), pch = c(21, 22), pt.bg = "white",
    col = "black", pt.cex = 1.2, bty = "n", horiz = TRUE, title = "Node type")
  if (length(sig_mods)) legend("topright", legend = paste0("M", sig_mods), pch = 15, pt.cex = 1.2,
    col = pal_sig, bty = "n", cex = 0.7, title = "ORA-significant modules")
  dev.off()

  list(g = g, edges_group = edges_group, top_proteins = top_proteins, top_traits = top_traits,
       prot_trait_comp = prot_trait_comp, trait_cluster_comp = trait_cluster_comp, fisher_df = fisher_df)
}
