#!/usr/bin/env Rscript
# =============================================================================
# make_example_data.R  --  write the public example datasets as plain files
# =============================================================================
# Creates input files in exactly the format a user would upload, from two
# classic Bioconductor datasets:
#
#   data/example/airway/
#     airway_counts.csv     63,677 genes (Ensembl IDs) x 8 samples
#     airway_metadata.csv   sample_id, condition (control / dexamethasone), cell_line
#     Human airway smooth muscle cells from 4 donors (cell lines), each treated
#     with dexamethasone or left untreated (Himes et al. 2014, PLoS One;
#     Bioconductor package 'airway'). Paired design: ~ cell_line + condition.
#
#   data/example/pasilla/
#     pasilla_counts.tsv    14,599 genes (FlyBase IDs) x 7 samples
#     pasilla_metadata.csv  sample_id, condition (untreated / treated), type
#     Drosophila S2 cells with pasilla knock-down (Brooks et al. 2011, Genome
#     Research; Bioconductor package 'pasilla'). 'type' (single-read vs
#     paired-end library) is a technical covariate: ~ type + condition.
#     The original annotation names samples "treated1fb" while the counts use
#     "treated1"; this script fixes that, which is exactly the kind of mismatch
#     step 01 is designed to catch.
#
#   data/example/gse60450/
#     gse60450_counts.tsv     27,179 genes (NCBI Entrez IDs) x 12 samples, plus the
#                             gene "Length" column, as published on GEO
#     gse60450_metadata.csv   sample_id, cell_type (basal / luminal),
#                             status (virgin / pregnant / lactating), lane
#     Mouse mammary gland basal and luminal cells from virgin, pregnant and
#     lactating mice, 2 replicates each (Fu et al. 2015, Nature Cell Biology;
#     GEO GSE60450). The dataset of the Bioconductor RNAseq123 workflow.
#     Kept close to the GEO file on purpose (Entrez IDs, extra Length column)
#     because that is what a biologist downloads. Sample names are shortened
#     to the library IDs (MCL1-DG ...). The counts file is downloaded from GEO,
#     or read from data/raw/GSE60450_Lactation-GenewiseCounts.txt.gz if it is there.
#     Note: lane is the same as cell type (basal = L002, luminal = L001), so
#     lane cannot be used as a covariate for cell type.
#
# Usage
#   Rscript scripts/make_example_data.R            # all datasets
#   Rscript scripts/make_example_data.R airway     # one dataset
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "scripts/x"), ".."))
})
which_sets <- commandArgs(trailingOnly = TRUE)
if (!length(which_sets)) which_sets <- c("airway", "pasilla", "gse60450")

if ("airway" %in% which_sets) {
  if (!requireNamespace("airway", quietly = TRUE)) {
    stop("Package 'airway' is not installed. Run: Rscript scripts/00_setup_environment.R")
  }
  suppressPackageStartupMessages(library(SummarizedExperiment))
  data("airway", package = "airway", envir = environment())
  out <- file.path(REPO_ROOT, "data", "example", "airway")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  cnt <- assay(airway, "counts")
  cd <- as.data.frame(colData(airway))
  meta <- data.frame(
    sample_id = rownames(cd),
    condition = ifelse(as.character(cd$dex) == "trt", "dexamethasone", "control"),
    cell_line = as.character(cd$cell),
    run = as.character(cd$Run),
    stringsAsFactors = FALSE
  )
  utils::write.csv(data.frame(gene_id = rownames(cnt), cnt, check.names = FALSE),
                   file.path(out, "airway_counts.csv"), row.names = FALSE, quote = FALSE)
  utils::write.csv(meta, file.path(out, "airway_metadata.csv"), row.names = FALSE, quote = FALSE)
  message("Wrote ", out, " (", nrow(cnt), " genes x ", ncol(cnt), " samples)")
}

if ("pasilla" %in% which_sets) {
  ext <- system.file("extdata", package = "pasilla")
  if (!nzchar(ext)) stop("Package 'pasilla' is not installed. Run: Rscript scripts/00_setup_environment.R")
  out <- file.path(REPO_ROOT, "data", "example", "pasilla")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  cnt <- utils::read.delim(file.path(ext, "pasilla_gene_counts.tsv"), check.names = FALSE)
  ann <- utils::read.csv(file.path(ext, "pasilla_sample_annotation.csv"))
  meta <- data.frame(
    sample_id = sub("fb$", "", ann$file),
    condition = as.character(ann$condition),
    type = as.character(ann$type),
    stringsAsFactors = FALSE
  )
  names(cnt)[1] <- "gene_id"
  utils::write.table(cnt, file.path(out, "pasilla_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.csv(meta, file.path(out, "pasilla_metadata.csv"), row.names = FALSE, quote = FALSE)
  message("Wrote ", out, " (", nrow(cnt), " genes x ", ncol(cnt) - 1, " samples)")
}

if ("gse60450" %in% which_sets) {
  raw <- file.path(REPO_ROOT, "data", "raw", "GSE60450_Lactation-GenewiseCounts.txt.gz")
  if (!file.exists(raw)) {
    dir.create(dirname(raw), recursive = TRUE, showWarnings = FALSE)
    url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/suppl/GSE60450_Lactation-GenewiseCounts.txt.gz"
    ok <- tryCatch(utils::download.file(url, raw, mode = "wb", quiet = TRUE) == 0, error = function(e) FALSE)
    if (!ok) {
      unlink(raw)
      stop("Could not download ", url, ". Download it in a browser and save it as ", raw)
    }
  }
  out <- file.path(REPO_ROOT, "data", "example", "gse60450")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  cnt <- utils::read.delim(gzfile(raw), check.names = FALSE)
  # "MCL1-DG_BC2CTUACXX_ACTTGA_L002_R1" -> sample "MCL1-DG", lane "L002"
  long <- names(cnt)[-(1:2)]
  ids <- sub("_.*$", "", long)
  names(cnt)[-(1:2)] <- ids
  code <- substr(ids, 6, 7)   # D = basal, L = luminal; G-L / A-F = mouse status
  meta <- data.frame(
    sample_id = ids,
    cell_type = ifelse(substr(code, 1, 1) == "D", "basal", "luminal"),
    status = c(G = "virgin", H = "virgin", I = "pregnant", J = "pregnant", K = "lactating", L = "lactating",
               A = "virgin", B = "virgin", C = "pregnant", D = "pregnant", E = "lactating", F = "lactating")[substr(code, 2, 2)],
    lane = sub("^.*_(L00[0-9])_.*$", "\\1", long),
    stringsAsFactors = FALSE
  )
  utils::write.table(cnt, file.path(out, "gse60450_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.csv(meta, file.path(out, "gse60450_metadata.csv"), row.names = FALSE, quote = FALSE)
  message("Wrote ", out, " (", nrow(cnt), " genes x ", length(ids), " samples)")
}
