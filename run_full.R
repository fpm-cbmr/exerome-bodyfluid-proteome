# =============================================================================
# run_full.R
#
# Full pipeline rebuild mode: analysis scripts plus figure/supplement rendering.
#
# Run:
#   Rscript run_full.R
# Optional flags:
#   RUN_ANALYSIS_14=0 Rscript run_full.R   # skip heavy cross-fluid network rebuild
#   RUN_ANALYSIS_15=1 Rscript run_full.R   # include OpenTargets coloc generation (analysis 15)
# =============================================================================
source(here::here("R/run_helpers.R"))

run_14 <- tolower(Sys.getenv("RUN_ANALYSIS_14", "1")) %in% c("1", "true", "yes")
run_15 <- tolower(Sys.getenv("RUN_ANALYSIS_15", "0")) %in% c("1", "true", "yes")

ensure_dirs(c(
  here::here("results", "period"),
  here::here("data", "network"),
  here::here("figures"),
  here::here("source_data")
))

# Core files required for full rebuild (00-13 and 15).
required_core <- c(
  here::here("data-raw", "olink", "exerome_olink.parquet"),
  here::here("data-raw", "body_fluid_proteomics", "urine.tsv"),
  here::here("data-raw", "body_fluid_proteomics", "saliva.tsv"),
  here::here("data-raw", "body_fluid_proteomics", "plasma.tsv"),
  here::here("data-raw", "body_fluid_proteomics", "metadata_urine_new.xlsx"),
  here::here("data-raw", "body_fluid_proteomics", "metadata_saliva.xlsx"),
  here::here("data-raw", "body_fluid_proteomics", "metadata_plasma_new.xlsx"),
  here::here("data-raw", "clinical_data.xlsx"),
  here::here("data-raw", "clinical_data", "glucose_cytokine_exercise_response.xlsx"),
  here::here("data-raw", "clinical_data", "covariates_msd.xlsx"),
  here::here("data-raw", "validation_exercise_mode", "validation_olink_Plate1.parquet"),
  here::here("data-raw", "validation_exercise_mode", "validation_olink_Plate2.parquet"),
  here::here("data-raw", "validation_exercise_mode", "validation_metadata.xlsx"),
  here::here("data-raw", "validation_clinical_dataa.xlsx"),
  here::here("data-raw", "pilot_clinical_dataa.xlsx"),
  here::here("data-raw", "hpa_24.tsv"),
  here::here("data-raw", "supplementary_table_2.xlsx"),
  here::here("data-raw", "supplementary_table_3.xlsx"),
  here::here("data-raw", "validation_exercise_mode", "exercise_order.xlsx"),
  here::here("data-raw", "phewas_pqtl", "newest_coloc_results_resubmission.xlsx"),
  here::here("data-raw", "phewas_pqtl", "final_pqtl_gwas_with_corrected_parents.csv"),
  here::here("data-raw", "approved_drug_targets.csv"),
  here::here("data-raw", "exerkine_list.xlsx"),
  here::here("multiomics_II", "figures", "network_diseases_secreted_testing_Combined_ORA.xlsx")
)
assert_files_exist(required_core, context = "full rebuild inputs")

if (run_15) {
  required_ot <- c(
    here::here("coloc_parquet"),
    here::here("crediblesets_parquet"),
    here::here("part-00000-1e42fdda-f476-428f-8a26-2c8b5c9a3b1e-c000.snappy.parquet"),
    here::here("all_sig_proteins.csv")
  )
  assert_files_exist(required_ot, context = "analysis 15 OpenTargets inputs")
}

analysis <- c(
  "R/analysis/00_data_processing.R",
  "R/analysis/01_metabolic_cytokine_lmm.R",
  "R/analysis/02_saliva_lmm_clusters.R",
  "R/analysis/02_saliva_tsne.R",
  "R/analysis/03_urine_lmm_clusters.R",
  "R/analysis/03_urine_tsne.R",
  "R/analysis/03b_plasma_lcms_lmm_clusters.R",
  "R/analysis/04_plasma_lmm_clusters.R",
  "R/analysis/04_plasma_tsne.R",
  "R/analysis/05_saliva_hpa_enrichment.R",
  "R/analysis/05b_hpa_disease_categories.R",
  "R/analysis/06_saliva_cluster_ora.R",
  "R/analysis/07_plasma_cluster_ora.R",
  "R/analysis/07b_urine_cluster_ora.R",
  "R/analysis/08_plasma_hpa_enrichment.R",
  "R/analysis/09_urine_hpa_enrichment.R",
  "R/analysis/10_bodyfluid_data_container.R",
  "R/analysis/10b_replication_olink_processing.R",
  "R/analysis/11_replication_period_lmm.R",
  "R/analysis/12_replication_sextime.R",
  "R/analysis/13_package_replication_period.R",
  "R/analysis/13b_replication_figure_objects.R",
  "R/analysis/16_pqtl_coloc_networks.R"
)

if (run_14) analysis <- c(analysis, "R/analysis/14_crossfluid_correlation_networks.R")
if (run_15) analysis <- c(analysis, "R/analysis/15_opentargets_coloc_generation.R")

run_many(analysis)

# --- render figures + supplementary from the analysis outputs above ----------
figures <- c(
  "R/figures/figure1.R", "R/figures/figure2.R", "R/figures/figure3.R",
  "R/figures/figure4.R", "R/figures/figure5.R", "R/figures/figure6.R",
  "R/figures/extended_data1.R", "R/figures/extended_data2.R",
  "R/figures/extended_data3.R", "R/figures/extended_data3b.R",
  "R/figures/extended_data4.R"
)
supplementary <- c(
  "R/supplementary/supplementary_figure_1_saliva_dynamics.R",
  "R/supplementary/supplementary_figure_2_saliva_tissue.R",
  "R/supplementary/supplementary_figure_3_urine_dynamics.R",
  "R/supplementary/supplementary_figure_4_urine_tissue.R",
  "R/supplementary/supplementary_figure_5_plasma_lcms.R",
  "R/supplementary/supplementary_figure_5e_olink_lcms_exercise_highlight.R",
  "R/supplementary/supplementary_figure_6a_exerkine_zscore_panel.R",
  "R/supplementary/supplementary_figure_6b_msd_olink_crossplatform_validation.R",
  "R/supplementary/supplementary_figure_7_plasma_dynamics.R",
  "R/supplementary/supplementary_figure_8_plasma_examples.R",
  "R/supplementary/supplementary_figure_9_plasma_tissue.R"
)
# source-data workbooks, rebuilt from the committed data/ analysis objects
source_data <- c(
  "R/source_data/figure1_source_data.R", "R/source_data/figure2_source_data.R",
  "R/source_data/figure3_source_data.R", "R/source_data/figure4_source_data.R",
  "R/source_data/figure5_source_data.R", "R/source_data/figure6_source_data.R",
  "R/source_data/extended_data1_source_data.R",
  "R/source_data/extended_data2_source_data.R", "R/source_data/extended_data3_source_data.R",
  "R/source_data/extended_data4_source_data.R",
  "R/source_data/supplementary_figure_1_source_data.R", "R/source_data/supplementary_figure_2_source_data.R",
  "R/source_data/supplementary_figure_3_source_data.R", "R/source_data/supplementary_figure_4_source_data.R",
  "R/source_data/supplementary_figure_5_source_data.R", "R/source_data/supplementary_figure_6_source_data.R",
  "R/source_data/supplementary_figure_7_source_data.R", "R/source_data/supplementary_figure_8_source_data.R",
  "R/source_data/supplementary_figure_9_source_data.R"
)
run_many(c(figures, supplementary, source_data))
