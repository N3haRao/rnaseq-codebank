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
# Usage
#   Rscript scripts/make_example_data.R            # both datasets
#   Rscript scripts/make_example_data.R airway     # one dataset
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "scripts/x"), ".."))
})
which_sets <- commandArgs(trailingOnly = TRUE)
if (!length(which_sets)) which_sets <- c("airway", "pasilla")

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
