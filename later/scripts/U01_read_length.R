#!/usr/bin/env Rscript
# =============================================================================
# U01_read_length.R  --  read length of FASTQ files (needed by rMATS)
# =============================================================================
# R version of read_length.py. Instead of the first read only, it looks at the
# first --n_reads reads of every file (trimmed files have mixed lengths) and
# reports the most common and the longest length.
#
# Usage
#   Rscript later/scripts/U01_read_length.R --fastq a_R1.fastq.gz,b_R1.fastq.gz --out_dir results/read_length
#   Rscript later/scripts/U01_read_length.R --sample_sheet samples.tsv --fastq_dir fastq/ --out_dir results/read_length
#
# Outputs (in --out_dir):
#   read_length.txt   one number: the read length to give rMATS (--readLength)
#   read_length.json  per file: reads checked, min / most common / max length;
#                     variable_length = TRUE means use rMATS --variable-read-length
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
  make_option("--fastq", type = "character", help = "Comma-separated FASTQ files (.fastq / .fastq.gz)"),
  make_option("--sample_sheet", type = "character", help = "Sample sheet (file_name, sample, read_orientation)"),
  make_option("--fastq_dir", type = "character", default = ".", help = "Folder with the FASTQ files named in the sample sheet"),
  make_option("--n_reads", type = "integer", default = 10000L, help = "Reads checked per file [default %default]"),
  make_option("--out_dir", type = "character", help = "Output folder (required)")
), description = "Read length of FASTQ files."))
if (is_blank(opt[["out_dir"]])) fail("--out_dir is required.")

files <- if (!is_blank(opt[["fastq"]])) {
  as_list_param(opt[["fastq"]])
} else if (!is_blank(opt[["sample_sheet"]])) {
  ss <- read_sample_sheet(opt[["sample_sheet"]])
  vapply(ss$file_name, function(f) if (file.exists(f)) f else file.path(opt[["fastq_dir"]], f), "")
} else fail("Give --fastq or --sample_sheet.")
missing <- files[!file.exists(files)]
if (length(missing)) fail("FASTQ file(s) not found: ", paste(missing, collapse = ", "))

read_lengths <- function(path, n) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, n = 4 * n, warn = FALSE)
  if (!length(lines)) fail("Empty FASTQ file: ", path)
  hdr <- lines[seq(1, length(lines), by = 4)]
  if (!all(startsWith(hdr, "@"))) fail("Not a FASTQ file (record headers must start with '@'): ", path)
  nchar(lines[seq(2, length(lines), by = 4)])
}

per_file <- lapply(files, function(f) {
  len <- read_lengths(f, opt[["n_reads"]])
  tab <- table(len)
  list(file = basename(f), reads_checked = length(len), min = min(len), max = max(len),
       most_common = as.integer(names(tab)[which.max(tab)]))
})
most_common <- vapply(per_file, `[[`, integer(1), "most_common")
maxes <- vapply(per_file, `[[`, integer(1), "max")
mins <- vapply(per_file, `[[`, integer(1), "min")
read_length <- as.integer(names(which.max(table(most_common))))
variable <- any(mins != maxes) || length(unique(most_common)) > 1

dir.create(opt[["out_dir"]], showWarnings = FALSE, recursive = TRUE)
writeLines(as.character(read_length), file.path(opt[["out_dir"]], "read_length.txt"))
write_json_file(list(read_length = read_length, longest = max(maxes), variable_length = variable,
                     advice = if (variable) "Reads differ in length (e.g. after trimming): run rMATS with --variable-read-length." else "All reads have the same length.",
                     files = per_file), file.path(opt[["out_dir"]], "read_length.json"))
log_info(sprintf("Read length %d (longest %d)%s; %d file(s) checked.", read_length, max(maxes),
                 if (variable) ", lengths vary" else "", length(files)))
