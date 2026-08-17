# =============================================================================
# Analysis 07b — Urine over-representation analysis (ORA), overall + per cluster.
#
# GO:BP / KEGG / Reactome over-representation with g:Profiler (gost), run right
# after the navmix clustering (analysis 03): once on all exercise-regulated urine
# proteins (overall), and once per temporal cluster. Queries the live g:Profiler
# database. (Supplementary Fig. 3f-h.)
#
# Input:  data/res.urine.linear.rda
# Output: data/res.enrich.urine.rda, data/res.enrich.urine.cluster.rda
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(stringr); library(purrr)
                   library(gprofiler2); library(data.table) })
load(here("data", "res.urine.linear.rda"))

gene_symbols <- function(x) unique(na.omit(unlist(str_split(x, "\\|"))))
ora <- function(genes, bg, srcs) gost(query = genes, organism = "hsapiens",
  ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
  measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05,
  correction_method = "fdr", domain_scope = "annotated", custom_bg = bg, numeric_ns = "",
  sources = srcs, as_short_link = FALSE)

prot.back <- gene_symbols(res.urine.linear$hgnc_symbol)

# ---- overall (KEGG + Reactome) ----------------------------------------------
res.enrich.urine <- ora(gene_symbols(res.urine.linear$hgnc_symbol[res.urine.linear$fdr.aov < 0.05]),
                        prot.back, c("KEGG", "REAC"))$result %>%
  mutate(fc = (intersection_size / term_size) / (query_size / effective_domain_size))
save(res.enrich.urine, file = here("data", "res.enrich.urine.rda"))

# ---- per temporal cluster (KEGG + Reactome + GO:BP) --------------------------
cl.label.exerome <- tibble(cluster = 1:5,
  label = c("Decrease, <baseline at 24h", "Decrease, >baseline at 24h",
            "Increase, <baseline at 24h", "Increase, >baseline at 24h", "Stochastic pattern"),
  cluster.paper = paste0("C", 1:5))
res.enrich.urine.cluster <- map_df(1:5, ~{
  g <- ora(gene_symbols(res.urine.linear$hgnc_symbol[res.urine.linear$fdr.aov < 0.05 &
                                                     res.urine.linear$cluster == .x]),
           prot.back, c("KEGG", "REAC", "GO:BP"))
  if (!is.null(g$result)) tibble(cluster = .x, result = g$result) else NULL
})
res.enrich.urine.cluster <- as.data.table(left_join(cl.label.exerome, res.enrich.urine.cluster, by = "cluster"))
res.enrich.urine.cluster[, fc := (result.intersection_size / result.term_size) /
                           (result.query_size / result.effective_domain_size)]
res.enrich.urine.cluster <- as.data.frame(res.enrich.urine.cluster)
save(res.enrich.urine.cluster, file = here("data", "res.enrich.urine.cluster.rda"))

cat(sprintf("Analysis 07b done: res.enrich.urine %d overall + res.enrich.urine.cluster %d cluster terms\n",
            nrow(res.enrich.urine), nrow(res.enrich.urine.cluster)))
