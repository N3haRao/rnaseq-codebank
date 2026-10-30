#!/usr/bin/env Rscript
# =============================================================================
# U03_tximport_counts.R  --  Salmon transcripts -> gene counts for the pipeline
# =============================================================================
# The bridge between U02 (Salmon) and the main pipeline, which starts from a
# gene counts table. Uses tximport (Soneson et al. 2015) with transcript-to-gene
# mapping from the GTF the Salmon index was built from.
#
# countsFromAbundance = "lengthScaledTPM" (default) gives counts corrected for
# transcript-length changes between samples, suitable as a plain counts table
# for DESeq2 (the code bank rounds them; round_counts: true). Use "no" to keep
# Salmon's estimated counts.
#
# Usage
#   Rscript later/scripts/U03_tximport_counts.R --salmon_dir results/salmon --gtf ref/genes.gtf \
#       --sample_sheet samples.tsv --out_dir results/counts
#
# Outputs (in --out_dir):
#   counts.csv           gene_id + one column per sample  -> the main pipeline's counts file
#   samples.csv          one row per sample (sample sheet without file columns) -> its sample table
#   gene_annotation.csv  gene_id, gene_name
#   tximport_summary.json
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "later/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "later", "R", "later_utils.R"))
suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--salmon_dir", type = "character", help = "U02 output folder (<sample>/quant.sf) (required)"),
  make_option("--gtf", type = "character", help = "GTF the Salmon index was built from"),
  make_option("--tx2gene", type = "character", help = "Instead of --gtf: CSV with transcript_id, gene_id[, gene_name]"),
  make_option("--sample_sheet", type = "character", help = "Sample sheet used for U02 (sample order and sample table)"),
  make_option("--counts_from_abundance", type = "character", default = "lengthScaledTPM",
              help = "no | scaledTPM | lengthScaledTPM [default %default]"),
  make_option("--ignore_tx_version", action = "store_true", default = FALSE, help = "Match ENST00000456328.2 to ENST00000456328"),
  make_option("--out_dir", type = "character", help = "Output folder (required)")
), description = "Gene counts from Salmon output with tximport."))
for (k in c("salmon_dir", "out_dir")) if (is_blank(opt[[k]])) fail("--", k, " is required.")
if (!requireNamespace("tximport", quietly = TRUE)) fail("Package 'tximport' is not installed (BiocManager::install('tximport')).")
if (!opt[["counts_from_abundance"]] %in% c("no", "scaledTPM", "lengthScaledTPM")) fail("--counts_from_abundance must be no, scaledTPM or lengthScaledTPM.")

# ---- Samples ---------------------------------------------------------------------------
quants <- list.files(opt[["salmon_dir"]], pattern = "^quant\\.sf$", recursive = TRUE, full.names = TRUE)
quants <- setNames(quants, basename(dirname(quants)))
sheet <- NULL
if (!is_blank(opt[["sample_sheet"]])) {
  ss <- read_sample_sheet(opt[["sample_sheet"]])
  sheet <- ss[!duplicated(ss$sample), setdiff(names(ss), c("file_name", "read_orientation")), drop = FALSE]
  missing <- setdiff(sheet$sample, names(quants))
  if (length(missing)) fail("No quant.sf for sample(s): ", paste(missing, collapse = ", "), " (did U02 finish for them?)")
  quants <- quants[sheet$sample]
}
if (!length(quants)) fail("No quant.sf files under ", opt[["salmon_dir"]])

# ---- Transcript -> gene ----------------------------------------------------------------------
tx2gene <- if (!is_blank(opt[["tx2gene"]])) {
  t <- read_table_file(opt[["tx2gene"]])
  if (!all(c("transcript_id", "gene_id") %in% names(t))) fail("--tx2gene needs columns transcript_id and gene_id.")
  t
} else if (!is_blank(opt[["gtf"]])) {
  tx2gene_from_gtf(opt[["gtf"]])
} else fail("Give --gtf or --tx2gene.")
log_info(sprintf("%d transcripts in %d genes.", nrow(tx2gene), length(unique(tx2gene$gene_id))))

# GENCODE transcript FASTA names look like ENST...|ENSG...|...; keep the first field
first_q <- data.table::fread(quants[[1]], data.table = FALSE, select = "Name")$Name
strip <- grepl("\\|", first_q[1])
if (opt[["ignore_tx_version"]]) tx2gene$transcript_id <- sub("\\.\\d+$", "", tx2gene$transcript_id)
tx_names <- if (strip) sub("\\|.*$", "", first_q) else first_q
if (opt[["ignore_tx_version"]]) tx_names <- sub("\\.\\d+$", "", tx_names)
matched <- mean(tx_names %in% tx2gene$transcript_id)
if (matched < 0.5) fail(sprintf("Only %.0f%% of Salmon's transcript names are in the GTF. Is it the GTF the index was built from? Try --ignore_tx_version.", 100 * matched))

txi <- tximport::tximport(quants, type = "salmon", tx2gene = tx2gene[, c("transcript_id", "gene_id")],
                          countsFromAbundance = opt[["counts_from_abundance"]], ignoreTxVersion = opt[["ignore_tx_version"]],
                          ignoreAfterBar = strip)
counts <- txi$counts
dir.create(opt[["out_dir"]], showWarnings = FALSE, recursive = TRUE)
write_table(data.frame(gene_id = rownames(counts), round(counts, 3), check.names = FALSE), file.path(opt[["out_dir"]], "counts.csv"))
ann <- tx2gene[!duplicated(tx2gene$gene_id), c("gene_id", "gene_name")]
write_table(ann[ann$gene_id %in% rownames(counts), ], file.path(opt[["out_dir"]], "gene_annotation.csv"))
if (!is.null(sheet)) write_table(sheet, file.path(opt[["out_dir"]], "samples.csv"))
write_json_file(list(n_genes = nrow(counts), n_samples = ncol(counts), samples = colnames(counts),
                     counts_from_abundance = opt[["counts_from_abundance"]], transcripts_matched_percent = round(100 * matched, 1),
                     tximport_version = as.character(utils::packageVersion("tximport")),
                     next_step = "Run the main pipeline with counts_file: counts.csv and metadata_file: samples.csv (round_counts: true)."),
                file.path(opt[["out_dir"]], "tximport_summary.json"))
log_info(sprintf("Gene counts: %d genes x %d samples -> %s", nrow(counts), ncol(counts), file.path(opt[["out_dir"]], "counts.csv")))
