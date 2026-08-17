# =============================================================================
# Analysis 07 — Plasma Olink over-representation analysis (ORA), overall + per cluster.
#
# GO:BP / KEGG / Reactome over-representation with g:Profiler (gost), run right
# after the navmix clustering (analysis 04): once on all exercise-regulated plasma
# proteins (overall), and once per temporal cluster. Queries the live g:Profiler
# database. (Supplementary Fig. 7g + Fig. 3 panel e.)
#
# Input:  data/res.olink.linear.rda
# Output: data/res.enrich.olink.rda, data/res.enrich.olink.cluster.rda
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(stringr); library(purrr)
                   library(gprofiler2); library(data.table) })
load(here("data", "res.olink.linear.rda"))

gene_symbols <- function(x) unique(na.omit(unlist(str_split(x, "\\|"))))
ora <- function(genes, bg, srcs) gost(query = genes, organism = "hsapiens",
  ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
  measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05,
  correction_method = "fdr", domain_scope = "annotated", custom_bg = bg, numeric_ns = "",
  sources = srcs, as_short_link = FALSE)

prot.back <- gene_symbols(res.olink.linear$hgnc_symbol)

# ---- overall (KEGG + Reactome) ----------------------------------------------
res.enrich.olink <- ora(gene_symbols(res.olink.linear$hgnc_symbol[res.olink.linear$fdr.aov < 0.05]),
                        prot.back, c("KEGG", "REAC"))$result %>%
  mutate(fc = (intersection_size / term_size) / (query_size / effective_domain_size))
save(res.enrich.olink, file = here("data", "res.enrich.olink.rda"))

# ---- per temporal cluster (KEGG + Reactome + GO:BP) --------------------------
cl.label.exerome <- tibble(cluster = 1:6,
  label = c("Decrease, over correction at 1h", "Increase, to baseline at 3h",
            "Increase, over correction", "Increase, slow decrease to 24h",
            "Increase, elevated at 24h", "Stochastic patterns"),
  cluster.paper = paste0("C", 1:6))
res.enrich.olink.cluster <- map_df(1:6, ~{
  g <- ora(gene_symbols(res.olink.linear$hgnc_symbol[res.olink.linear$fdr.aov < 0.05 &
                                                     res.olink.linear$cluster == .x]),
           prot.back, c("KEGG", "REAC", "GO:BP"))
  if (!is.null(g$result)) tibble(cluster = .x, result = g$result) else NULL
})
res.enrich.olink.cluster <- as.data.table(left_join(cl.label.exerome, res.enrich.olink.cluster, by = "cluster"))
res.enrich.olink.cluster[, fc := (result.intersection_size / result.term_size) /
                           (result.query_size / result.effective_domain_size)]
res.enrich.olink.cluster <- as.data.frame(res.enrich.olink.cluster)
save(res.enrich.olink.cluster, file = here("data", "res.enrich.olink.cluster.rda"))

cat(sprintf("Analysis 07 done: res.enrich.olink %d overall + res.enrich.olink.cluster %d cluster terms\n",
            nrow(res.enrich.olink), nrow(res.enrich.olink.cluster)))
