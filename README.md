# Systemic proteomic dynamics of acute exercise reveal genetic links to health and disease

This repository contains code that was used to analyze and produce all outputs from the published paper.
However, it requires dedicated input files, which are not sharable due to GDPR, e.g., individual level proteomics and clinical data. 

## Organization

- `R/analysis/`: analysis pipeline scripts.
- `R/figures/`: main and extended-data figure scripts.
- `R/supplementary/`: supplementary figure scripts.
- `R/functions_loading.R`, `R/package_loading.R`, `R/figure_defaults.R`: shared utilities/config.

Directories such as `data/`, `data-raw/`, `results/`, `benchmark/`, `release/`,
are not present. 

## Analysis Script Map

The analysis pipeline is organized as numbered scripts in `R/analysis/`.

1. `R/analysis/00_data_processing.R`
Builds per-fluid NPX/protein tables from raw files.

2. `R/analysis/01_metabolic_cytokine_lmm.R`
Fits LMMs for metabolic markers/cytokines in the discovery data.

3. `R/analysis/02_saliva_lmm_clusters.R`
Runs saliva LMMs and temporal clustering for exercise-responsive proteins.

4. `R/analysis/02_saliva_tsne.R`
Computes saliva sample t-SNE coordinates.

5. `R/analysis/03_urine_lmm_clusters.R`
Runs urine LMMs and temporal clustering for exercise-responsive proteins.

6. `R/analysis/03_urine_tsne.R`
Computes urine sample t-SNE coordinates.

7. `R/analysis/03b_plasma_lcms_lmm_clusters.R`
Runs plasma LC-MS LMMs and temporal clustering; writes `res.plasma.linear` and the figure-ready plasma exerome/protein-label objects.

8. `R/analysis/04_plasma_lmm_clusters.R`
Runs plasma Olink LMMs, temporal clustering, and variance decomposition.

9. `R/analysis/04_plasma_tsne.R`
Computes plasma sample t-SNE coordinates.

10. `R/analysis/05_saliva_hpa_enrichment.R`
Performs saliva cluster tissue/cell-type enrichment using HPA annotations.

11. `R/analysis/05b_hpa_disease_categories.R`
Builds `data_hpa_categorized` from HPA disease annotations plus approved drug-target labels.

12. `R/analysis/06_saliva_cluster_ora.R`
Runs live g:Profiler ORA for saliva (overall + per cluster) and writes both ORA objects.

13. `R/analysis/07_plasma_cluster_ora.R`
Runs live g:Profiler ORA for plasma Olink (overall + per cluster) and writes both ORA objects.

14. `R/analysis/07b_urine_cluster_ora.R`
Runs live g:Profiler ORA for urine (overall + per cluster) and writes both ORA objects.

15. `R/analysis/08_plasma_hpa_enrichment.R`
Performs plasma cluster tissue/cell-type enrichment using HPA annotations.

16. `R/analysis/09_urine_hpa_enrichment.R`
Performs urine cluster tissue/cell-type enrichment using HPA annotations.

17. `R/analysis/10_bodyfluid_data_container.R`
Builds combined cross-fluid container (`se_bodyfluid`) for network/correlation analyses.

18. `R/analysis/10b_replication_olink_processing.R`
Processes replication Olink NPX/protein-label/normalized objects for downstream replication models.

19. `R/analysis/11_replication_period_lmm.R`
Fits replication-cohort LMMs adjusted for exercise-period randomization.

20. `R/analysis/12_replication_sextime.R`
Runs corrected replication three-way (`sex x time x mode`) sensitivity models.

21. `R/analysis/13_package_replication_period.R`
Packages replication outputs into figure-ready objects in `data/`.

22. `R/analysis/13b_replication_figure_objects.R`
Builds Fig. 5 / Extended Data 2-3 replication figure objects (matches, labels, enrichments).

23. `R/analysis/14_crossfluid_correlation_networks.R`
Builds cross-fluid correlation/network objects for discovery analyses.

24. `R/analysis/15_opentargets_coloc_generation.R`
Optional OpenTargets colocalization generation workflow (external/HPC step; output committed).

25. `R/analysis/16_pqtl_coloc_networks.R`
Builds pQTL colocalization integration/network artifacts from the committed coloc output.


## Contact

For questions about missing restricted inputs or how to run the analysis, please reach out: nigel.kurgan@sund.ku.dk
