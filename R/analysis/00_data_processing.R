# =============================================================================
# Analysis 00 — Raw data processing: vendor exports -> per-fluid NPX matrices.
#
# Turns the raw proteomics exports into the wide per-fluid abundance tables and
# protein-annotation tables that the rest of the pipeline (02/03/04 onward)
# starts from. 
#
# Raw inputs (place under data-raw/):
#   data-raw/olink/exerome_olink.parquet            Olink Explore NPX export (plasma)
#   data-raw/body_fluid_proteomics/plasma.tsv       LC-MS (Spectronaut) report, plasma
#   data-raw/body_fluid_proteomics/saliva.tsv       LC-MS (Spectronaut) report, saliva
#   data-raw/body_fluid_proteomics/urine.tsv        LC-MS (Spectronaut) report, urine
#   data-raw/body_fluid_proteomics/metadata_*.xlsx  per-fluid sample sheets
#   data-raw/clinical_data.xlsx                     subject covariates
#
# Outputs (data/):
#   olink_npx.data  / olink.prot.label              plasma Olink NPX + annotation
#   saliva_npx_data / prot.label.saliva             saliva LC-MS NPX + annotation
#   urine_npx_data  / prot.label.urine              urine  LC-MS NPX + annotation
#   plasma_npx_data / prot.label.plasma             plasma LC-MS NPX + annotation
#   (+ se_*_new / *_processed intermediates)
# =============================================================================
suppressMessages({
  library(here); library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(readr); library(readxl); library(arrow)
  library(AnnotationDbi); library(org.Hs.eg.db); library(biomaRt)
  library(PhosR); library(SummarizedExperiment)
})

TIME_FROM_LABEL <- function(id) dplyr::case_when(
  grepl("PRE",  id) ~ -1, grepl("T0",  id) ~ 0,  grepl("T30", id) ~ 0.5,
  grepl("T60",  id) ~ 1,  grepl("T180", id) ~ 3, grepl("T24h", id) ~ 24,
  TRUE ~ NA_real_)

# map gene symbols -> Ensembl gene id + UniProt (org.Hs.eg.db / BioMart), then
# enrich with gene coordinates, HGNC symbol, description and exon ids. This yields
# the 12-column prot.label.* annotation that every downstream script consumes
# (the same BioMart attribute set used in the original per-fluid analyses).
annotate_symbols <- function(symbols) {
  lab <- data.frame(protein = symbols, Assay = symbols, stringsAsFactors = FALSE)
  lab$ensembl_ids <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = lab$Assay,
                                           column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
  mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ids  <- unique(stats::na.omit(lab$ensembl_ids))
  # one UniProt / external-gene-name row per gene
  map  <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "uniprotswissprot"),
                         filters = "ensembl_gene_id", values = ids, mart = mart)
  lab  <- merge(lab, map, by.x = "ensembl_ids", by.y = "ensembl_gene_id", all.x = TRUE) %>%
    dplyr::distinct(Assay, .keep_all = TRUE)
  # gene coordinates + HGNC symbol + description + exon ids (exon-expanded, as published)
  annotation <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "chromosome_name", "strand", "hgnc_symbol",
                   "start_position", "end_position", "description", "ensembl_exon_id"),
    filters = "ensembl_gene_id", values = ids, mart = mart)
  merge(lab, annotation, by.x = "ensembl_ids", by.y = "ensembl_gene_id", all.x = TRUE)
}

# read one LC-MS Spectronaut report and produce a log2, median-scaled SummarizedExperiment
process_lcms <- function(tsv, sample_cols, first_col, metadata_xlsx, drop_col = "PG.ProteinGroups") {
  raw <- readr::read_tsv(here(tsv), show_col_types = FALSE)
  new_names <- sample_cols
  names(raw)[first_col:(first_col + length(new_names) - 1)] <- new_names

  proc <- raw %>%
    dplyr::mutate(PG.Genes = ifelse(PG.Genes == "", PG.ProteinGroups, PG.Genes)) %>%
    dplyr::mutate_all(~ ifelse(is.nan(.), NA, .)) %>%
    dplyr::mutate_at(first_col:(first_col + length(new_names) - 1), as.numeric) %>%
    dplyr::mutate(PG.Genes = make.names(PG.Genes, unique = TRUE)) %>%
    tibble::column_to_rownames("PG.Genes") %>%
    dplyr::select(-dplyr::any_of(c("PG.ProteinGroups", "PG.UniProtIds", "PG.MolecularWeight"))) %>%
    log2() %>%
    tibble::rownames_to_column("entry") %>%
    dplyr::mutate(entry = sub("\\..*", "", entry), entry = make.names(entry, unique = TRUE)) %>%
    tibble::column_to_rownames("entry")

  metadata <- readxl::read_excel(here(metadata_xlsx))
  # `remove` flags samples to drop; some per-fluid sheets omit it when nothing is dropped
  if (!"remove" %in% names(metadata)) metadata$remove <- NA_character_
  drop <- metadata %>% dplyr::filter(remove == "x") %>% dplyr::pull(sample_id)
  filt <- proc %>% dplyr::select(-dplyr::any_of(drop))
  scaled <- filt %>% PhosR::selectOverallPercent(percent = 0.5) %>% PhosR::medianScaling()
  meta_keep <- metadata %>% dplyr::filter(is.na(remove) | remove != "x") %>% dplyr::select(-remove)
  se <- PhosR::PhosphoExperiment(assay = scaled, colData = meta_keep)
  list(processed = proc, se = se, metadata = metadata)
}

# turn a per-fluid SummarizedExperiment into a wide npx table (rows = samples)
se_to_npx <- function(se, metadata, drop_extra = character()) {
  data <- as.data.frame(SummarizedExperiment::assay(se))
  tdat <- as.data.frame(t(data)); tdat$sample_id <- rownames(tdat)
  tdat <- dplyr::left_join(tdat, metadata, by = "sample_id") %>%
    dplyr::filter(is.na(remove) | remove != "x")
  tdat %>%
    dplyr::select(sample_id, time_point = time_numeric, dplyr::everything(),
                  -dplyr::any_of(c("remove", "sample", "sample_id", "plate",
                                   "plate_position", "time", "biological_id", drop_extra))) %>%
    dplyr::mutate(label = as.factor(subject), time = as.factor(time_point)) %>%
    dplyr::select(-subject, -time_point)
}

npx_prot_label <- function(npx) {
  cols <- setdiff(names(npx), c("time", "label"))
  annotate_symbols(cols)
}

# ---- Plasma Olink ------------------------------------------------------------
df <- arrow::read_parquet(here("data-raw/olink/exerome_olink.parquet"), as_data_frame = TRUE) %>%
  dplyr::mutate(subject = stringr::str_extract(SampleID, "^EX\\d+"),
                time_label = TIME_FROM_LABEL(SampleID)) %>%
  dplyr::filter(!is.na(subject), AssayType == "assay")

olink_npx.data <- tidyr::pivot_wider(df, id_cols = c(time_label, subject),
                                     names_from = OlinkID, values_from = NPX) %>%
  dplyr::mutate(label = as.factor(subject), time = as.factor(time_label)) %>%
  dplyr::select(-subject, -time_label) %>%
  dplyr::mutate(label = stringr::str_remove(label, "^EX"))
save(olink_npx.data, file = here("data", "olink_npx.data.rda"))

olink.prot.label <- df %>% dplyr::distinct(OlinkID, Assay, UniProt)
olink.prot.label$ensembl_ids <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = olink.prot.label$Assay,
                                                      column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
save(olink.prot.label, file = here("data", "olink.prot.label.rda"))

# ---- Urine LC-MS (75 samples: S_1..S_75, data columns start at 5) -----------
u <- process_lcms("data-raw/body_fluid_proteomics/urine.tsv", paste0("S_", 1:75), 5,
                  "data-raw/body_fluid_proteomics/metadata_urine_new.xlsx")
urine_processed <- u$processed; se_urine_new <- u$se
save(urine_processed, file = here("data", "urine_processed.rda"))
save(se_urine_new,    file = here("data", "se_urine_new.rda"))
urine_npx_data  <- se_to_npx(u$se, u$metadata)
save(urine_npx_data, file = here("data", "urine_npx_data.rda"))
prot.label.urine <- npx_prot_label(urine_npx_data)
save(prot.label.urine, file = here("data", "prot.label.urine.rda"))

# ---- Saliva LC-MS (S_213..S_267, S_269..S_349; data columns start at 3) -----
s <- process_lcms("data-raw/body_fluid_proteomics/saliva.tsv",
                  c(paste0("S_", 213:267), paste0("S_", 269:349)), 3,
                  "data-raw/body_fluid_proteomics/metadata_saliva.xlsx")
saliva_processed <- s$processed; se_saliva_new <- s$se
save(saliva_processed, file = here("data", "saliva_processed.rda"))
save(se_saliva_new,    file = here("data", "se_saliva_new.rda"))
saliva_npx_data  <- se_to_npx(s$se, s$metadata)
save(saliva_npx_data, file = here("data", "saliva_npx_data.rda"))
prot.label.saliva <- npx_prot_label(saliva_npx_data)
save(prot.label.saliva, file = here("data", "prot.label.saliva.rda"))

# ---- Plasma LC-MS (S_76..S_212; data columns start at 3) --------------------
p <- process_lcms("data-raw/body_fluid_proteomics/plasma.tsv", paste0("S_", 76:212), 3,
                  "data-raw/body_fluid_proteomics/metadata_plasma_new.xlsx")
plasma_processed <- p$processed; se_plasma_new <- p$se
save(plasma_processed, file = here("data", "plasma_processed.rda"))
save(se_plasma_new,    file = here("data", "se_plasma_new.rda"))
plasma_npx_data  <- se_to_npx(p$se, p$metadata)
save(plasma_npx_data, file = here("data", "plasma_npx_data.rda"))
prot.label.plasma <- npx_prot_label(plasma_npx_data)
save(prot.label.plasma, file = here("data", "prot.label.plasma.rda"))

cat(sprintf("Analysis 00 done: olink %d, urine %d, saliva %d, plasma %d proteins.\n",
            ncol(olink_npx.data) - 2, nrow(urine_processed),
            nrow(saliva_processed), nrow(plasma_processed)))
