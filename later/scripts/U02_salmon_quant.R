#!/usr/bin/env Rscript
# =============================================================================
# U02_salmon_quant.R  --  transcript quantification with Salmon
# =============================================================================
# R version of salmon_alignment.py. Changes: single-end samples are supported
# (a sample with only an R1 / SE file), every sample's exit code is checked,
# and the library type Salmon detects (-l A) is turned into the strandedness
# that rMATS needs (strandedness.txt), so no separate strand-detection step.
#
# Sample sheet (TSV or CSV), one row per FASTQ file:
#   file_name             sample  condition  read_orientation  batch
#   CT1_R1.fastq.gz       CT1     CT         R1                1
#   CT1_R2.fastq.gz       CT1     CT         R2                1
#
# Usage
#   Rscript later/scripts/U02_salmon_quant.R --sample_sheet samples.tsv --fastq_dir fastq/ \
#       --index ref/salmon_index --out_dir results/salmon [--threads 8 --jobs 2 --bootstraps 30]
#
# Outputs (in --out_dir):
#   <sample>/quant.sf ...      Salmon's own output per sample
#   salmon_run_summary.csv     status, seconds, mapping rate, detected library type per sample
#   salmon_run_summary.json    the same plus parameters and Salmon version
#   strandedness.txt           unstranded | forward | reverse (agreed by all samples)
#   logs/<sample>.log          Salmon's console output
# Exit code 1 if any sample failed.
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
  make_option("--sample_sheet", type = "character", help = "Sample sheet (required)"),
  make_option("--fastq_dir", type = "character", default = ".", help = "Folder with the FASTQ files"),
  make_option("--index", type = "character", help = "Salmon index folder (required)"),
  make_option("--out_dir", type = "character", help = "Output folder (required)"),
  make_option("--sample_col", type = "character", default = "sample"),
  make_option("--file_col", type = "character", default = "file_name"),
  make_option("--orientation_col", type = "character", default = "read_orientation"),
  make_option("--lib_type", type = "character", default = "A", help = "Salmon library type; A = detect [default %default]"),
  make_option("--threads", type = "integer", default = 8L, help = "Threads per sample [default %default]"),
  make_option("--jobs", type = "integer", default = 1L, help = "Samples run at the same time [default %default]"),
  make_option("--bootstraps", type = "integer", default = 30L, help = "Bootstrap replicates [default %default]"),
  make_option("--no_gc_bias", action = "store_true", default = FALSE),
  make_option("--no_seq_bias", action = "store_true", default = FALSE),
  make_option("--no_pos_bias", action = "store_true", default = FALSE),
  make_option("--extra", type = "character", help = "Extra Salmon arguments, space-separated"),
  make_option("--timeout", type = "integer", default = 7200L, help = "Seconds per sample [default %default]"),
  make_option("--salmon", type = "character", help = "Path to salmon (default: on PATH)"),
  make_option("--min_mapping_rate", type = "double", default = 50, help = "Warn below this mapping rate in percent [default %default]"),
  make_option("--dry_run", action = "store_true", default = FALSE, help = "Print the commands only")
), description = "Salmon quantification of FASTQ files."))
for (k in c("sample_sheet", "index", "out_dir")) if (is_blank(opt[[k]])) fail("--", k, " is required.")

salmon <- need_tool("salmon", opt[["salmon"]])
if (!dir.exists(opt[["index"]]) || !file.exists(file.path(opt[["index"]], "info.json"))) {
  fail("The Salmon index '", opt[["index"]], "' is missing or incomplete (no info.json). Build it with: salmon index -t transcripts.fa -i <index>")
}
ss <- read_sample_sheet(opt[["sample_sheet"]], opt[["sample_col"]], opt[["file_col"]], opt[["orientation_col"]])
samples <- pair_fastqs(ss, opt[["fastq_dir"]], opt[["sample_col"]], opt[["file_col"]], opt[["orientation_col"]])
dir.create(opt[["out_dir"]], showWarnings = FALSE, recursive = TRUE)
out_dir <- normalizePath(opt[["out_dir"]])
version <- tool_version(salmon)
log_info(sprintf("%s; %d samples (%d paired-end, %d single-end).", version, nrow(samples),
                 sum(!is.na(samples$fastq_2)), sum(is.na(samples$fastq_2))))

run_one <- function(i) {
  s <- samples[i, ]
  reads <- if (is.na(s$fastq_2)) c("-r", shQuote(s$fastq_1)) else c("-1", shQuote(s$fastq_1), "-2", shQuote(s$fastq_2))
  args <- c("quant", "-i", shQuote(normalizePath(opt[["index"]])), "-l", opt[["lib_type"]], reads,
            "-o", shQuote(file.path(out_dir, s$sample)), "--validateMappings",
            "--numBootstraps", opt[["bootstraps"]], "--threads", opt[["threads"]],
            if (!opt[["no_gc_bias"]] && !is.na(s$fastq_2)) "--gcBias",   # GC bias needs paired-end reads
            if (!opt[["no_seq_bias"]]) "--seqBias", if (!opt[["no_pos_bias"]]) "--posBias",
            if (!is_blank(opt[["extra"]])) strsplit(opt[["extra"]], " +")[[1]])
  t0 <- Sys.time()
  code <- run_tool(salmon, args, file.path(out_dir, "logs", paste0(s$sample, ".log")), opt[["timeout"]], opt[["dry_run"]])
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  quant <- file.path(out_dir, s$sample, "quant.sf")
  ok <- code == 0 && (opt[["dry_run"]] || file.exists(quant))
  meta <- read_json_file(file.path(out_dir, s$sample, "aux_info", "meta_info.json"))
  lib <- read_json_file(file.path(out_dir, s$sample, "lib_format_counts.json"))
  list(sample = s$sample, status = if (ok) "success" else "failed", exit_code = code, seconds = secs,
       read_type = if (is.na(s$fastq_2)) "single" else "paired",
       mapping_rate = if (!is.null(meta$percent_mapped)) round(meta$percent_mapped, 2) else NA,
       n_reads = meta$num_processed %||% NA, library_type = lib$expected_format %||% NA,
       error = if (ok) "" else if (code != 0) sprintf("salmon exit code %d (see logs/%s.log)", code, s$sample) else "quant.sf not written")
}
results <- if (opt[["jobs"]] > 1 && .Platform$OS.type == "unix") {
  parallel::mclapply(seq_len(nrow(samples)), run_one, mc.cores = min(opt[["jobs"]], nrow(samples)))
} else lapply(seq_len(nrow(samples)), run_one)
res <- do.call(rbind, lapply(results, as.data.frame))
write_table(res, file.path(out_dir, "salmon_run_summary.csv"))

# ---- Library type -> strandedness --------------------------------------------------------
lt <- unique(na.omit(res$library_type))
strand <- unique(unname(SALMON_STRAND[lt]))
strand_note <- if (length(strand) == 1 && !is.na(strand)) {
  writeLines(strand, file.path(out_dir, "strandedness.txt"))
  sprintf("All samples: library type %s = %s.", paste(lt, collapse = "/"), strand)
} else if (!length(lt)) "No library type recorded (dry run, or Salmon failed)." else
  sprintf("Samples disagree on the library type (%s); strandedness.txt not written. Check the libraries.", paste(lt, collapse = ", "))
log_info(strand_note)
low <- res$sample[!is.na(res$mapping_rate) & res$mapping_rate < opt[["min_mapping_rate"]]]
if (length(low)) log_warn(sprintf("Low mapping rate (< %g%%): %s. Check the index species and adapter trimming.", opt[["min_mapping_rate"]], paste(low, collapse = ", ")))

write_json_file(list(salmon_version = version, index = normalizePath(opt[["index"]]), parameters = opt[setdiff(names(opt), "help")],
                     strandedness = if (length(strand) == 1) strand else NULL, strandedness_note = strand_note,
                     low_mapping_rate = as.list(low), samples = results),
                file.path(out_dir, "salmon_run_summary.json"))
n_fail <- sum(res$status == "failed")
if (n_fail) fail(sprintf("%d of %d samples failed: %s", n_fail, nrow(res), paste(res$sample[res$status == "failed"], collapse = ", ")))
log_info(sprintf("Salmon finished for all %d samples -> %s", nrow(res), out_dir))
