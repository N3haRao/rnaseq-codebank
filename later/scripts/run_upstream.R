#!/usr/bin/env Rscript
# =============================================================================
# run_upstream.R  --  from FASTQ files to everything downstream, in one command
# =============================================================================
# Runs the later/ scripts in order with one config file, stops at the first
# failure, and writes upstream_summary.json with every output path. Meant for
# the Globus upload integration: FASTQ files + a sample sheet in, counts and
# splicing / isoform results out.
#
#   qc           U00 fastp trimming and QC        -> trimmed FASTQ + samples_trimmed.tsv
#   read_length  U01                               -> read_length.txt
#   star         U06 STAR alignment (+ index)      -> BAM files, STAR gene counts, strandedness
#   salmon       U02 Salmon (+ index if missing)   -> quant.sf per sample
#   gene_counts  U03 tximport                      -> counts.csv + samples.csv for the main pipeline
#   rmats        U04 rMATS, each group vs control  -> rMATS tables + splicing configs
#   splicing     splicing module per comparison
#   sashimi      U05 per comparison
#   isoforms     isoform module (DEI + DUI) per comparison
#
# Usage
#   Rscript later/scripts/run_upstream.R --config upstream.yaml
#
# Config (YAML):
#   sample_sheet: samples.tsv          # file_name, sample, read_orientation, <condition column>, ...
#   fastq_dir: fastq/
#   out_dir: results/
#   condition_column: genotype
#   control: WT
#   genome_fasta: ref/genome.fa        # for the STAR index (built once)
#   gtf: ref/genes.gtf
#   star_index: ref/star_index
#   salmon_index: ref/salmon_index     # built from transcripts_fasta when missing
#   transcripts_fasta: ref/transcripts.fa
#   threads: 8
#   steps: [qc, read_length, star, salmon, gene_counts, rmats, splicing, sashimi, isoforms]
#   splicing_settings: {}              # extra settings for the splicing module, e.g. {delta_psi_cutoff: 0.2}
#   isoform_settings: {}               # extra settings for the isoform module
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
  make_option("--config", type = "character", help = "Upstream config YAML (required)")
), description = "FASTQ to counts, splicing and isoform results."))
if (is_blank(opt[["config"]]) || !file.exists(opt[["config"]])) fail("--config is required (see the header of this script).")
cfg <- yaml::read_yaml(opt[["config"]])
ALL_STEPS <- c("qc", "read_length", "star", "salmon", "gene_counts", "rmats", "splicing", "sashimi", "isoforms")
steps <- as_list_param(cfg$steps %||% ALL_STEPS)
bad <- setdiff(steps, ALL_STEPS)
if (length(bad)) fail("Unknown step(s): ", paste(bad, collapse = ", "), ". Steps: ", paste(ALL_STEPS, collapse = ", "))
for (k in c("sample_sheet", "out_dir")) if (is_blank(cfg[[k]])) fail("The config needs '", k, "'.")
abs_path <- function(p) if (is_blank(p)) NULL else normalizePath(p, mustWork = FALSE)
out <- abs_path(cfg$out_dir); dir.create(out, recursive = TRUE, showWarnings = FALSE)
threads <- as.character(cfg$threads %||% 4)
sheet <- abs_path(cfg$sample_sheet); fastq_dir <- abs_path(cfg$fastq_dir %||% ".")
rscript <- file.path(R.home("bin"), "Rscript")
summary <- list(started = timestamp_now(), config = cfg, steps = list())

step <- function(name, script, args) {
  log_info(sprintf("== %s ==", name))
  log_file <- file.path(out, "logs", paste0(name, ".log"))
  dir.create(dirname(log_file), showWarnings = FALSE)
  t0 <- Sys.time()
  code <- system2(rscript, c(shQuote(file.path(REPO_ROOT, script)), args), stdout = log_file, stderr = log_file)
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  summary$steps[[name]] <<- list(status = if (code == 0) "success" else "failed", seconds = secs, log = log_file)
  write_json_file(summary, file.path(out, "upstream_summary.json"))
  if (code != 0) {
    cat(paste0("    | ", tail(readLines(log_file, warn = FALSE), 12), "\n"), sep = "")
    fail(sprintf("Step '%s' failed (see %s).", name, log_file))
  }
  log_info(sprintf("   done in %.0f s", secs))
}
newest_run <- function(dir, before) setdiff(list.dirs(dir, recursive = FALSE), before)[1]
# extra module settings from the config -> --set "key=value;key2=[a, b]"
set_args <- function(extra) {
  if (!length(extra)) return(character(0))
  val <- function(v) if (length(v) > 1 || is.list(v)) paste0("[", paste(unlist(v), collapse = ", "), "]") else as.character(v)
  c("--set", shQuote(paste(sprintf("%s=%s", names(extra), vapply(extra, val, "")), collapse = ";")))
}

if ("qc" %in% steps) {
  step("qc", "later/scripts/U00_fastq_qc_trim.R", c("--sample_sheet", shQuote(sheet), "--fastq_dir", shQuote(fastq_dir),
                                                    "--out_dir", shQuote(file.path(out, "trimmed")), "--threads", threads))
  sheet <- file.path(out, "trimmed", "samples_trimmed.tsv"); fastq_dir <- file.path(out, "trimmed")
}
if ("read_length" %in% steps) {
  step("read_length", "later/scripts/U01_read_length.R", c("--sample_sheet", shQuote(sheet), "--fastq_dir", shQuote(fastq_dir),
                                                           "--out_dir", shQuote(file.path(out, "read_length"))))
}
read_length <- if (file.exists(file.path(out, "read_length", "read_length.txt"))) file.path(out, "read_length", "read_length.txt") else as.character(cfg$read_length %||% 101)
if ("star" %in% steps) {
  step("star", "later/scripts/U06_star_align.R", c("--sample_sheet", shQuote(sheet), "--fastq_dir", shQuote(fastq_dir),
                                                   "--star_index", shQuote(abs_path(cfg$star_index %||% file.path(out, "star_index"))),
                                                   if (!is_blank(cfg$genome_fasta)) c("--genome_fasta", shQuote(abs_path(cfg$genome_fasta))),
                                                   if (!is_blank(cfg$gtf)) c("--gtf", shQuote(abs_path(cfg$gtf))),
                                                   "--read_length", shQuote(read_length), "--threads", threads,
                                                   "--out_dir", shQuote(file.path(out, "star"))))
}
if ("salmon" %in% steps) {
  idx <- abs_path(cfg$salmon_index %||% file.path(out, "salmon_index"))
  if (!file.exists(file.path(idx, "info.json"))) {
    if (is_blank(cfg$transcripts_fasta)) fail("No Salmon index at ", idx, " and no transcripts_fasta to build one.")
    log_info("== salmon_index ==")
    salmon <- need_tool("salmon")
    if (run_tool(salmon, c("index", "-t", shQuote(abs_path(cfg$transcripts_fasta)), "-i", shQuote(idx), "-p", threads),
                 file.path(out, "logs", "salmon_index.log")) != 0) fail("Building the Salmon index failed (logs/salmon_index.log).")
  }
  step("salmon", "later/scripts/U02_salmon_quant.R", c("--sample_sheet", shQuote(sheet), "--fastq_dir", shQuote(fastq_dir),
                                                       "--index", shQuote(idx), "--out_dir", shQuote(file.path(out, "salmon")),
                                                       "--threads", threads, "--bootstraps", as.character(cfg$bootstraps %||% 30)))
}
if ("gene_counts" %in% steps) {
  if (is_blank(cfg$gtf)) fail("gene_counts needs 'gtf' in the config.")
  step("gene_counts", "later/scripts/U03_tximport_counts.R", c("--salmon_dir", shQuote(file.path(out, "salmon")), "--gtf", shQuote(abs_path(cfg$gtf)),
                                                               "--sample_sheet", shQuote(sheet), "--out_dir", shQuote(file.path(out, "counts"))))
}
# strandedness: STAR's call, otherwise Salmon's, otherwise the config
strand_file <- c(file.path(out, "star", "strandedness.txt"), file.path(out, "salmon", "strandedness.txt"))
strand <- if (any(file.exists(strand_file))) strand_file[file.exists(strand_file)][1] else (cfg$strandedness %||% "unstranded")
paired <- any(toupper(read_sample_sheet(sheet)$read_orientation) %in% c("R2", "2"))
comparisons <- character(0)
if ("rmats" %in% steps) {
  for (k in c("condition_column", "control", "gtf")) if (is_blank(cfg[[k]])) fail("rmats needs '", k, "' in the config.")
  step("rmats", "later/scripts/U04_rmats.R", c("--sample_sheet", shQuote(sheet), "--condition_column", cfg$condition_column,
                                               "--control", shQuote(cfg$control), "--bam_dir", shQuote(file.path(out, "star", "bams")),
                                               "--gtf", shQuote(abs_path(cfg$gtf)), "--read_length", shQuote(read_length),
                                               "--strandedness", shQuote(strand), "--read_type", if (paired) "paired" else "single",
                                               "--threads", threads, "--out_dir", shQuote(file.path(out, "rmats"))))
}
if (dir.exists(file.path(out, "rmats"))) comparisons <- list.dirs(file.path(out, "rmats"), recursive = FALSE)
runs_dir <- file.path(out, "runs"); dir.create(runs_dir, showWarnings = FALSE)
summary$splicing_runs <- list()
for (cmp in comparisons) {
  name <- basename(cmp)
  if ("splicing" %in% steps) {
    set <- set_args(cfg$splicing_settings)
    before <- list.dirs(runs_dir, recursive = FALSE)
    step(paste0("splicing_", name), "splicing/scripts/run_splicing.R", c("--config", shQuote(file.path(cmp, "splicing_config.yaml")),
                                                                         "--runs_dir", shQuote(runs_dir), set))
    srun <- newest_run(runs_dir, before)
    summary$splicing_runs[[name]] <- srun
    if ("sashimi" %in% steps) {
      step(paste0("sashimi_", name), "later/scripts/U05_sashimi_plots.R", c("--splicing_run", shQuote(srun), "--comparison_dir", shQuote(cmp),
                                                                            "--gtf", shQuote(abs_path(cfg$gtf))))
    }
  }
}
if ("isoforms" %in% steps) {
  for (k in c("condition_column", "control", "gtf")) if (is_blank(cfg[[k]])) fail("isoforms needs '", k, "' in the config.")
  smp <- file.path(out, "counts", "samples.csv")
  if (!file.exists(smp)) fail("isoforms needs the gene_counts step (it writes the sample table).")
  groups <- setdiff(unique(read.csv(smp, check.names = FALSE)[[cfg$condition_column]]), cfg$control)
  summary$isoform_runs <- list()
  set <- set_args(cfg$isoform_settings)
  for (g in groups) {
    before <- list.dirs(runs_dir, recursive = FALSE)
    step(paste0("isoforms_", make.names(g)), "isoforms/scripts/run_isoforms.R",
         c("--quant_dir", shQuote(file.path(out, "salmon")), "--metadata_file", shQuote(smp), "--gtf_file", shQuote(abs_path(cfg$gtf)),
           "--condition", cfg$condition_column, "--reference_level", shQuote(cfg$control), "--test_level", shQuote(g),
           "--run_name", make.names(paste0(cfg$control, "_vs_", g, "_isoforms")), "--runs_dir", shQuote(runs_dir), set))
    summary$isoform_runs[[g]] <- newest_run(runs_dir, before)
  }
}
summary$finished <- timestamp_now()
summary$main_pipeline_inputs <- list(
  counts_salmon = if (file.exists(file.path(out, "counts", "counts.csv"))) file.path(out, "counts", "counts.csv") else NULL,
  counts_star = if (file.exists(file.path(out, "star", "counts_star.csv"))) file.path(out, "star", "counts_star.csv") else NULL,
  sample_table = if (file.exists(file.path(out, "counts", "samples.csv"))) file.path(out, "counts", "samples.csv") else
    if (file.exists(file.path(out, "star", "samples.csv"))) file.path(out, "star", "samples.csv") else NULL)
write_json_file(summary, file.path(out, "upstream_summary.json"))
log_info("All steps finished. Summary: ", file.path(out, "upstream_summary.json"))
