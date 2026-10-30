#!/usr/bin/env Rscript
# =============================================================================
# U00_fastq_qc_trim.R  --  read quality control and adapter / quality trimming
# =============================================================================
# The first step for raw FASTQ files (missing from the lab's scripts). Runs fastp
# (Chen et al. 2018) on every sample: adapter detection and removal, trimming of
# low-quality ends, removal of reads that become too short, and a QC report.
# Paired-end and single-end samples are both supported.
#
# Usage
#   Rscript later/scripts/U00_fastq_qc_trim.R --sample_sheet samples.tsv --fastq_dir fastq/ --out_dir results/trimmed
#
# Outputs (in --out_dir):
#   <sample>_R1.trimmed.fastq.gz (+ _R2)   trimmed reads
#   reports/<sample>.html / .json           fastp's report per sample
#   qc_summary.csv / .json                  reads before / after, % kept, Q30, adapters, duplication, flags
#   samples_trimmed.tsv                     the sample sheet pointing at the trimmed files (input for U01-U06)
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
  make_option("--sample_sheet", type = "character", help = "Sample sheet (file_name, sample, read_orientation, ...) (required)"),
  make_option("--fastq_dir", type = "character", default = ".", help = "Folder with the FASTQ files"),
  make_option("--out_dir", type = "character", help = "Output folder (required)"),
  make_option("--threads", type = "integer", default = 4L, help = "Threads per sample [default %default]"),
  make_option("--jobs", type = "integer", default = 1L, help = "Samples at the same time [default %default]"),
  make_option("--min_length", type = "integer", default = 25L, help = "Drop reads shorter than this after trimming [default %default]"),
  make_option("--qualified_quality", type = "integer", default = 15L, help = "Base quality counted as good [default %default]"),
  make_option("--cut_tail_quality", type = "integer", default = 20L, help = "Trim 3' ends below this mean quality (sliding window) [default %default]"),
  make_option("--adapter_r1", type = "character", help = "Adapter sequence for read 1 (default: detected)"),
  make_option("--adapter_r2", type = "character", help = "Adapter sequence for read 2 (default: detected)"),
  make_option("--min_kept_percent", type = "double", default = 80, help = "Flag samples keeping fewer reads [default %default]"),
  make_option("--min_q30_percent", type = "double", default = 75, help = "Flag samples with fewer Q30 bases after trimming [default %default]"),
  make_option("--fastp", type = "character", help = "Path to fastp (default: on PATH)"),
  make_option("--dry_run", action = "store_true", default = FALSE)
), description = "FASTQ quality control and trimming with fastp."))
for (k in c("sample_sheet", "out_dir")) if (is_blank(opt[[k]])) fail("--", k, " is required.")
fastp <- need_tool("fastp", opt[["fastp"]])
ss <- read_sample_sheet(opt[["sample_sheet"]])
samples <- pair_fastqs(ss, opt[["fastq_dir"]])
dir.create(file.path(opt[["out_dir"]], "reports"), recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(opt[["out_dir"]])
version <- tool_version(fastp)
log_info(sprintf("%s; %d samples.", version, nrow(samples)))

run_one <- function(i) {
  s <- samples[i, ]
  pe <- !is.na(s$fastq_2)
  o1 <- file.path(out_dir, paste0(s$sample, "_R1.trimmed.fastq.gz"))
  o2 <- if (pe) file.path(out_dir, paste0(s$sample, "_R2.trimmed.fastq.gz")) else NA
  rep_json <- file.path(out_dir, "reports", paste0(s$sample, ".json"))
  args <- c("-i", shQuote(s$fastq_1), "-o", shQuote(o1),
            if (pe) c("-I", shQuote(s$fastq_2), "-O", shQuote(o2), "--detect_adapter_for_pe"),
            if (!is_blank(opt[["adapter_r1"]])) c("--adapter_sequence", opt[["adapter_r1"]]),
            if (pe && !is_blank(opt[["adapter_r2"]])) c("--adapter_sequence_r2", opt[["adapter_r2"]]),
            "--qualified_quality_phred", opt[["qualified_quality"]], "--cut_tail", "--cut_tail_mean_quality", opt[["cut_tail_quality"]],
            "--length_required", opt[["min_length"]], "--thread", opt[["threads"]],
            "--json", shQuote(rep_json), "--html", shQuote(file.path(out_dir, "reports", paste0(s$sample, ".html"))),
            "--report_title", shQuote(paste("fastp report:", s$sample)))
  code <- run_tool(fastp, args, file.path(out_dir, "logs", paste0(s$sample, ".log")), dry_run = opt[["dry_run"]])
  ok <- code == 0 && (opt[["dry_run"]] || file.exists(o1))
  j <- read_json_file(rep_json)
  before <- j$summary$before_filtering$total_reads %||% NA
  after <- j$summary$after_filtering$total_reads %||% NA
  list(sample = s$sample, status = if (ok) "success" else "failed", exit_code = code, read_type = if (pe) "paired" else "single",
       reads_before = before, reads_after = after,
       percent_kept = if (is.numeric(before) && before > 0) round(100 * after / before, 2) else NA,
       q30_percent_after = if (!is.null(j)) round(100 * j$summary$after_filtering$q30_rate, 2) else NA,
       mean_length_after = j$summary$after_filtering$read1_mean_length %||% NA,
       adapter_trimmed_reads = j$adapter_cutting$adapter_trimmed_reads %||% NA,
       duplication_percent = if (!is.null(j$duplication$rate)) round(100 * j$duplication$rate, 2) else NA,
       trimmed_1 = o1, trimmed_2 = o2)
}
res_list <- if (opt[["jobs"]] > 1 && .Platform$OS.type == "unix") parallel::mclapply(seq_len(nrow(samples)), run_one, mc.cores = min(opt[["jobs"]], nrow(samples))) else
  lapply(seq_len(nrow(samples)), run_one)
res <- do.call(rbind, lapply(res_list, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
res$flags <- trimws(paste(ifelse(!is.na(res$percent_kept) & res$percent_kept < opt[["min_kept_percent"]], sprintf("kept < %g%% of reads", opt[["min_kept_percent"]]), ""),
                          ifelse(!is.na(res$q30_percent_after) & res$q30_percent_after < opt[["min_q30_percent"]], sprintf("Q30 < %g%%", opt[["min_q30_percent"]]), "")))
write_table(res[, setdiff(names(res), c("trimmed_1", "trimmed_2"))], file.path(out_dir, "qc_summary.csv"))

# sample sheet for the next steps, pointing at the trimmed files
extra <- samples[, setdiff(names(samples), c("sample", "fastq_1", "fastq_2")), drop = FALSE]
trimmed <- do.call(rbind, lapply(seq_len(nrow(res)), function(i) {
  r1 <- data.frame(file_name = basename(res$trimmed_1[i]), sample = res$sample[i], read_orientation = "R1", extra[i, , drop = FALSE], check.names = FALSE)
  if (is.na(res$trimmed_2[i])) r1 else rbind(r1, transform(r1, file_name = basename(res$trimmed_2[i]), read_orientation = "R2"))
}))
write.table(trimmed, file.path(out_dir, "samples_trimmed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
flagged <- res$sample[nzchar(res$flags)]
if (length(flagged)) log_warn("Check these samples: ", paste(sprintf("%s (%s)", flagged, res$flags[nzchar(res$flags)]), collapse = "; "))
write_json_file(list(fastp_version = version, parameters = opt[setdiff(names(opt), "help")], flagged = as.list(flagged),
                     next_step = "Use samples_trimmed.tsv with --fastq_dir <this folder> for U01, U02 and U06.",
                     samples = res_list), file.path(out_dir, "qc_summary.json"))
n_fail <- sum(res$status == "failed")
if (n_fail) fail(sprintf("fastp failed for %d sample(s): %s", n_fail, paste(res$sample[res$status == "failed"], collapse = ", ")))
log_info(sprintf("Trimmed %d samples; median %.1f%% of reads kept -> %s", nrow(res), median(res$percent_kept, na.rm = TRUE), out_dir))
