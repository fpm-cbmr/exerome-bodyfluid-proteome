# =============================================================================
# Analysis 06 — Saliva over-representation analysis (ORA), overall + per cluster.
#
# GO:BP / KEGG / Reactome over-representation with g:Profiler (gost), run right
# after the navmix clustering (analysis 02): once on all exercise-regulated saliva
# proteins (overall), and once per temporal cluster. Queries the live g:Profiler
# database. (Supplementary Fig. 1h + Fig. 2 panel d.)
#
# Input:  data/res.saliva.linear.rda
# Output: data/res.enrich.saliva.rda, data/res.enrich.saliva.cluster.rda
# =============================================================================
suppressMessages({ library(here); library(dplyr); library(stringr); library(purrr)
                   library(gprofiler2); library(data.table) })
load(here("data", "res.saliva.linear.rda"))

gene_symbols <- function(x) unique(na.omit(unlist(str_split(x, "\\|"))))
ora <- function(genes, bg, srcs) gost(query = genes, organism = "hsapiens",
  ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
  measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05,
  correction_method = "fdr", domain_scope = "annotated", custom_bg = bg, numeric_ns = "",
  sources = srcs, as_short_link = FALSE)

prot.back <- gene_symbols(res.saliva.linear$hgnc_symbol)   # background = all measured proteins

# ---- overall: all exercise-regulated proteins (KEGG + Reactome) --------------
res.enrich.saliva <- ora(gene_symbols(res.saliva.linear$hgnc_symbol[res.saliva.linear$fdr.aov < 0.05]),
                         prot.back, c("KEGG", "REAC"))$result %>%
  mutate(fc = (intersection_size / term_size) / (query_size / effective_domain_size))
save(res.enrich.saliva, file = here("data", "res.enrich.saliva.rda"))

# ---- per temporal cluster (KEGG + Reactome + GO:BP) --------------------------
cl.label.exerome <- tibble(cluster = 1:7,
  label = c("Increase, elevated at 24h", "Decrease, baseline at 24h",
            "Increase at 3h, elevated at 24h", "Decrease, reduced at 24h",
            "Increase, to baseline at 1-3h", "Sequential and sustained increase",
            "Stochastic patterns"),
  cluster.paper = paste0("C", 1:7))
res.enrich.saliva.cluster <- map_df(1:7, ~{
  g <- ora(gene_symbols(res.saliva.linear$hgnc_symbol[res.saliva.linear$fdr.aov < 0.05 &
                                                      res.saliva.linear$cluster == .x]),
           prot.back, c("KEGG", "REAC", "GO:BP"))
  if (!is.null(g$result)) tibble(cluster = .x, result = g$result) else NULL
})
res.enrich.saliva.cluster <- as.data.table(left_join(cl.label.exerome, res.enrich.saliva.cluster, by = "cluster"))
res.enrich.saliva.cluster[, fc := (result.intersection_size / result.term_size) /
                            (result.query_size / result.effective_domain_size)]
res.enrich.saliva.cluster <- as.data.frame(res.enrich.saliva.cluster)
save(res.enrich.saliva.cluster, file = here("data", "res.enrich.saliva.cluster.rda"))

cat(sprintf("Analysis 06 done: res.enrich.saliva %d overall + res.enrich.saliva.cluster %d cluster terms\n",
            nrow(res.enrich.saliva), nrow(res.enrich.saliva.cluster)))
