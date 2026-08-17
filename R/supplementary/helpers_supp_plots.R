# =============================================================================
# Shared plotting helpers for the supplementary figures
# =============================================================================
suppressMessages({ library(dplyr); library(ggplot2); library(ggpubr); library(ggrepel); library(writexl) })

# ---- source-data writer -----------------------------------------------------
# Writes/updates source_data/SourceData_SupplementaryFigure<fig>.xlsx. Merges with
# any sheets already in the file, so figures split across two scripts (S5, S6)
# both contribute to a single workbook regardless of run order.
write_supp_source_data <- function(sheets, fig) {
  f <- here::here("source_data", sprintf("SourceData_SupplementaryFigure%s.xlsx", fig))
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
  existing <- list()
  if (file.exists(f)) {
    sn <- readxl::excel_sheets(f)
    existing <- setNames(lapply(sn, function(s) as.data.frame(readxl::read_xlsx(f, sheet = s))), sn)
  }
  existing[names(sheets)] <- lapply(sheets, as.data.frame)
  writexl::write_xlsx(existing, f)
}

# ---- Disease-enrichment plot (incident-disease HR>1) ------------------------
# Signed -log10(P) on y, with point colour by cluster (grey if non-significant),
# size by odds ratio, and top diseases per cluster labelled.
# `df` must contain Disease, P_value, Odds_Ratio, Cluster, and padj.
determine_disease_order <- function(data, order_cluster = "5") {
  ord <- data %>% dplyr::filter(Cluster == order_cluster) %>%
    dplyr::arrange(P_value) %>% dplyr::pull(Disease) %>% unique()
  data %>% dplyr::mutate(Disease = factor(Disease, levels = ord))
}
prepare_plot_data <- function(data, cluster_color_palette) {
  data %>% dplyr::mutate(
    Log10P_signed = ifelse(Odds_Ratio > 1, -log10(P_value), log10(P_value)),
    Point_Color   = ifelse(padj < 0.05, as.character(Cluster), "Non-significant"),
    Label_Color   = ifelse(padj < 0.05, cluster_color_palette[as.character(Cluster)], "grey80"))
}
get_top2_per_cluster <- function(data) {
  data %>% dplyr::filter(padj < 0.05) %>% dplyr::group_by(Cluster) %>%
    dplyr::arrange(P_value) %>% dplyr::slice_head(n = 2) %>% dplyr::ungroup()
}
plot_disease_enrichment <- function(df, cluster_color_palette, title, order_cluster = "5") {
  df <- df %>% mutate(Cluster = as.character(Cluster)) %>%
    determine_disease_order(order_cluster) %>% prepare_plot_data(cluster_color_palette)
  top_data <- get_top2_per_cluster(df) %>% filter(Cluster != "1")
  ggplot(df, aes(x = Disease, y = Log10P_signed)) +
    geom_point(aes(color = Point_Color, size = Odds_Ratio), alpha = 0.7) +
    scale_color_manual(values = c(cluster_color_palette, "Non-significant" = "grey90"),
                       name = "Cluster", breaks = c(names(cluster_color_palette), "Non-significant")) +
    scale_size_continuous(range = c(0.5, 3), name = "Odds Ratio") +
    geom_text_repel(data = top_data, aes(x = Disease, y = Log10P_signed, label = Disease, color = Cluster),
                    inherit.aes = FALSE, size = 2, force = 40, box.padding = 0.3, max.overlaps = 50,
                    show.legend = FALSE) +
    labs(title = title, x = "Disease", y = "Signed -log10(P-value)") +
    theme_pubr() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "right", legend.key.size = unit(2, "mm"), text = element_text(size = 6))
}

# ---- Secreted-protein enrichment per cluster (Fisher) -----------------------
# Fisher enrichment of annotated secreted proteins per cluster.
# `ttp` is top_tissue_per_protein_* with cluster and Secretome.function columns.
plot_secretome_enrichment <- function(ttp, cluster_color_palette) {
  annotated_data <- ttp %>%
    dplyr::mutate(cluster.y = cluster,
                  annotated = ifelse(!is.na(Secretome.function) & Secretome.function != "",
                                     "Annotated", "Not Annotated"))
  results <- list()
  for (cl in unique(annotated_data$cluster.y)) {
    inc <- annotated_data %>% dplyr::filter(cluster.y == cl) %>% dplyr::count(annotated) %>%
      tidyr::complete(annotated = c("Annotated", "Not Annotated"), fill = list(n = 0)) %>%
      dplyr::rename(in_cluster_count = n)
    out <- annotated_data %>% dplyr::filter(cluster.y != cl) %>% dplyr::count(annotated) %>%
      tidyr::complete(annotated = c("Annotated", "Not Annotated"), fill = list(n = 0)) %>%
      dplyr::rename(not_in_cluster_count = n)
    ct <- inc %>% dplyr::left_join(out, by = "annotated") %>% dplyr::mutate(
      in_cluster_annotated = ifelse(annotated == "Annotated", in_cluster_count, 0),
      in_cluster_not_annotated = ifelse(annotated == "Not Annotated", in_cluster_count, 0),
      not_in_cluster_annotated = ifelse(annotated == "Annotated", not_in_cluster_count, 0),
      not_in_cluster_not_annotated = ifelse(annotated == "Not Annotated", not_in_cluster_count, 0)) %>%
      dplyr::select(in_cluster_annotated, in_cluster_not_annotated,
                    not_in_cluster_annotated, not_in_cluster_not_annotated) %>%
      dplyr::summarise_all(sum)
    m <- matrix(c(ct$in_cluster_annotated, ct$in_cluster_not_annotated,
                  ct$not_in_cluster_annotated, ct$not_in_cluster_not_annotated), nrow = 2, byrow = TRUE)
    ft <- stats::fisher.test(m)
    results[[length(results) + 1]] <- tibble::tibble(
      cluster.y = cl, p_value = ft$p.value, odds_ratio = unname(ft$estimate),
      in_cluster_annotated = ct$in_cluster_annotated)
  }
  results_df <- dplyr::bind_rows(results) %>%
    dplyr::mutate(p_adj = stats::p.adjust(p_value, "BH"), log_adj_p = -log10(p_adj)) %>%
    dplyr::filter(!is.na(cluster.y))
  ggplot(results_df, aes(x = odds_ratio, y = log_adj_p)) +
    geom_point(aes(fill = factor(cluster.y), size = in_cluster_annotated),
               shape = 21, color = "black", stroke = ifelse(results_df$p_adj < 0.05, 1.5, 0.5)) +
    scale_fill_manual(values = cluster_color_palette, name = NULL) +
    scale_size_continuous(name = "Secretome counts", range = c(1, 3), breaks = c(10, 50, 150)) +
    labs(title = NULL, x = "Odds Ratio", y = expression(-log[10](Adj.~p~val))) +
    theme_pubr() +
    theme(text = element_text(size = 6), legend.position = c(0.99, 0.7),
          legend.justification = c("right", "bottom"), legend.key.size = unit(2, "mm"),
          legend.background = element_blank(), legend.title = element_text(size = 6),
          legend.text = element_text(size = 6), plot.title = element_text(hjust = 0.5)) +
    guides(fill = "none", size = guide_legend(override.aes = list(stroke = 0.5)))
}
