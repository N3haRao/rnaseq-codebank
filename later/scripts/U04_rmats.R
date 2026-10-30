#!/usr/bin/env Rscript
# =============================================================================
# U04_rmats.R  --  run rMATS: each group against the control group
# =============================================================================
# R version of rmats_analysis.py (the rMATS part; figures and statistics
# afterwards are done by the splicing module, sashimi plots by U05). Changes:
#   * rMATS's exit code and output files are checked; a failed comparison stops
#     the script with exit code 1 instead of being skipped silently
#   * strandedness maps to valid rMATS library types (fr-unstranded,
#     fr-firststrand, fr-secondstrand), read from U02's strandedness.txt or given
#   * variable read lengths (from U01) switch on --variable-read-length
#   * each comparison writes splicing_config.yaml, ready for the splicing module:
#       Rscript splicing/scripts/run_splicing.R --config <out_dir>/<comparison>/splicing_config.yaml
#
# Usage
#   Rscript later/scripts/U04_rmats.R --sample_sheet samples.tsv --bam_dir bams/ \
#       --condition_column genotype --control WT --gtf ref/genes.gtf \
#       --read_length results/read_length/read_length.txt --strandedness results/salmon/strandedness.txt \
#       --read_type paired --out_dir results/rmats
#
# Outputs (in --out_dir/<control>_vs_<group>/):
#   tables/          rMATS output (SE.MATS.JCEC.txt, ...)
#   b1.txt, b2.txt   BAM lists (b1 = the group, b2 = the control)
#   splicing_config.yaml, rmats.log
# and rmats_summary.json in --out_dir.
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
  make_option("--sample_sheet", type = "character", help = "Samples with a condition column (required)"),
  make_option("--sample_col", type = "character", default = "sample"),
  make_option("--condition_column", type = "character", help = "Column with the groups (required)"),
  make_option("--control", type = "character", help = "Control group (required)"),
  make_option("--groups", type = "character", help = "Groups to compare with the control (default: all others)"),
  make_option("--bam_dir", type = "character", help = "Folder with <sample>.bam (or .sorted.bam, ...) (required)"),
  make_option("--gtf", type = "character", help = "Annotation GTF (required)"),
  make_option("--read_length", type = "character", help = "Read length, or U01's read_length.txt (required)"),
  make_option("--strandedness", type = "character", default = "unstranded",
              help = "unstranded | forward | reverse, or U02's strandedness.txt [default %default]"),
  make_option("--read_type", type = "character", help = "paired | single (required)"),
  make_option("--threads", type = "integer", default = 8L),
  make_option("--variable_read_length", action = "store_true", default = FALSE),
  make_option("--allow_clipping", action = "store_true", default = FALSE),
  make_option("--novel_ss", action = "store_true", default = FALSE, help = "Detect unannotated splice sites"),
  make_option("--cstat", type = "double", default = 0.0001, help = "rMATS splicing-difference cutoff of the test [default %default]"),
  make_option("--anchor_length", type = "integer", help = "Minimum junction anchor (rMATS default 1)"),
  make_option("--individual_counts", action = "store_true", default = FALSE),
  make_option("--rmats", type = "character", help = "Path to rmats.py (default: on PATH)"),
  make_option("--out_dir", type = "character", help = "Output folder (required)"),
  make_option("--dry_run", action = "store_true", default = FALSE)
), description = "rMATS for each group against the control."))
for (k in c("sample_sheet", "condition_column", "control", "bam_dir", "gtf", "read_length", "read_type", "out_dir")) {
  if (is_blank(opt[[k]])) fail("--", k, " is required.")
}
if (!opt[["read_type"]] %in% c("paired", "single")) fail("--read_type must be paired or single.")
rmats <- need_tool("rmats.py", opt[["rmats"]])
if (!file.exists(opt[["gtf"]])) fail("GTF not found: ", opt[["gtf"]])

# ---- Read length and strandedness ------------------------------------------------------
read_len <- if (grepl("^[0-9]+$", opt[["read_length"]])) as.integer(opt[["read_length"]]) else {
  if (!file.exists(opt[["read_length"]])) fail("Read length file not found: ", opt[["read_length"]])
  rl_json <- read_json_file(file.path(dirname(opt[["read_length"]]), "read_length.json"))
  if (isTRUE(rl_json$variable_length) && !opt[["variable_read_length"]]) {
    log_info("U01 found variable read lengths: using --variable-read-length.")
    opt[["variable_read_length"]] <- TRUE
  }
  as.integer(trimws(readLines(opt[["read_length"]], n = 1)))
}
if (is.na(read_len) || read_len < 20) fail("Read length must be a number of at least 20; got '", opt[["read_length"]], "'.")
strand <- if (opt[["strandedness"]] %in% names(RMATS_LIBTYPE)) opt[["strandedness"]] else {
  if (!file.exists(opt[["strandedness"]])) fail("--strandedness must be unstranded, forward, reverse or a file; got '", opt[["strandedness"]], "'.")
  trimws(readLines(opt[["strandedness"]], n = 1))
}
if (!strand %in% names(RMATS_LIBTYPE)) fail("Unknown strandedness '", strand, "' (use unstranded, forward or reverse).")
lib_type <- RMATS_LIBTYPE[[strand]]

# ---- Samples and BAM files ---------------------------------------------------------------------
ss <- read_table_file(opt[["sample_sheet"]])
for (k in c(opt[["sample_col"]], opt[["condition_column"]])) if (!k %in% names(ss)) fail("The sample sheet has no column '", k, "'.")
ss <- ss[!duplicated(ss[[opt[["sample_col"]]]]), ]
find_bam <- function(s) {
  cand <- file.path(opt[["bam_dir"]], paste0(s, c(".bam", ".sorted.bam", "_sorted.bam", ".aligned.bam", "_aligned.bam",
                                            "Aligned.sortedByCoord.out.bam", ".Aligned.sortedByCoord.out.bam")))
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(normalizePath(hit[1]))
  g <- Sys.glob(file.path(opt[["bam_dir"]], paste0("*", s, "*.bam")))
  if (length(g) == 1) return(normalizePath(g))
  if (length(g) > 1) fail(sprintf("Several BAM files match sample '%s': %s. Rename them to <sample>.bam.", s, paste(basename(g), collapse = ", ")))
  NA_character_
}
ss$bam <- vapply(ss[[opt[["sample_col"]]]], find_bam, "")
if (any(is.na(ss$bam))) fail("No BAM file for sample(s): ", paste(ss[[opt[["sample_col"]]]][is.na(ss$bam)], collapse = ", "), " in ", opt[["bam_dir"]])
groups <- unique(as.character(ss[[opt[["condition_column"]]]]))
if (!opt[["control"]] %in% groups) fail(sprintf("Control '%s' is not in column '%s' (groups: %s).", opt[["control"]], opt[["condition_column"]], paste(groups, collapse = ", ")))
tests <- if (!is_blank(opt[["groups"]])) as_list_param(opt[["groups"]]) else setdiff(groups, opt[["control"]])
bad <- setdiff(tests, groups)
if (length(bad)) fail("Unknown group(s): ", paste(bad, collapse = ", "))
version <- tool_version(rmats)
log_info(sprintf("%s; read length %d (%s reads), library %s; %d comparison(s) against %s.",
                 version, read_len, opt[["read_type"]], lib_type, length(tests), opt[["control"]]))

# ---- One rMATS run per comparison -----------------------------------------------------------------
dir.create(opt[["out_dir"]], showWarnings = FALSE, recursive = TRUE)
out_dir <- normalizePath(opt[["out_dir"]])
records <- list()
for (g in tests) {
  name <- make.names(paste0(opt[["control"]], "_vs_", g))
  cdir <- file.path(out_dir, name)
  tables <- file.path(cdir, "tables")
  dir.create(file.path(cdir, "tmp"), recursive = TRUE, showWarnings = FALSE)
  dir.create(tables, showWarnings = FALSE)
  test_s <- ss[ss[[opt[["condition_column"]]]] == g, ]
  ctrl_s <- ss[ss[[opt[["condition_column"]]]] == opt[["control"]], ]
  writeLines(paste(test_s$bam, collapse = ","), file.path(cdir, "b1.txt"))
  writeLines(paste(ctrl_s$bam, collapse = ","), file.path(cdir, "b2.txt"))
  args <- c("--b1", shQuote(file.path(cdir, "b1.txt")), "--b2", shQuote(file.path(cdir, "b2.txt")),
            "--gtf", shQuote(normalizePath(opt[["gtf"]])), "-t", opt[["read_type"]], "--readLength", read_len,
            "--libType", lib_type, "--nthread", opt[["threads"]], "--od", shQuote(tables), "--tmp", shQuote(file.path(cdir, "tmp")),
            "--cstat", opt[["cstat"]],
            if (opt[["variable_read_length"]]) "--variable-read-length", if (opt[["allow_clipping"]]) "--allow-clipping",
            if (opt[["novel_ss"]]) "--novelSS", if (opt[["individual_counts"]]) "--individual-counts",
            if (!is.null(opt[["anchor_length"]])) c("--anchorLength", opt[["anchor_length"]]))
  log_info(sprintf("%s: %d %s vs %d %s samples ...", name, nrow(test_s), g, nrow(ctrl_s), opt[["control"]]))
  code <- run_tool(rmats, args, file.path(cdir, "rmats.log"), dry_run = opt[["dry_run"]])
  produced <- file.exists(file.path(tables, "SE.MATS.JCEC.txt"))
  ok <- code == 0 && (opt[["dry_run"]] || produced)
  if (ok) unlink(file.path(cdir, "tmp"), recursive = TRUE)
  write_yaml_file(list(run_name = paste0(name, "_splicing"), rmats_dir = tables, count_type = "JCEC",
                       group_1_name = g, group_2_name = opt[["control"]],
                       group_1_samples = as.list(test_s[[opt[["sample_col"]]]]), group_2_samples = as.list(ctrl_s[[opt[["sample_col"]]]]),
                       rmats_version = sub("^.*v", "v", version)),
                  file.path(cdir, "splicing_config.yaml"),
                  header = "Config for the splicing module: Rscript splicing/scripts/run_splicing.R --config <this file>")
  records[[name]] <- list(comparison = name, group_1 = g, group_2 = opt[["control"]], n_1 = nrow(test_s), n_2 = nrow(ctrl_s),
                          status = if (ok) "success" else "failed", exit_code = code, tables = tables,
                          splicing_config = file.path(cdir, "splicing_config.yaml"))
  if (!ok) log_msg(sprintf("%s failed: %s", name, if (code != 0) "rMATS exit code " else "no SE.MATS.JCEC.txt written; see rmats.log"), code, level = "ERROR")
}
write_json_file(list(rmats_version = version, read_length = read_len, read_type = opt[["read_type"]], library_type = lib_type,
                     strandedness = strand, variable_read_length = opt[["variable_read_length"]], gtf = normalizePath(opt[["gtf"]]),
                     comparisons = unname(records)), file.path(out_dir, "rmats_summary.json"))
failed <- names(records)[vapply(records, function(r) r$status == "failed", logical(1))]
if (length(failed)) fail("rMATS failed for: ", paste(failed, collapse = ", "), ". See <comparison>/rmats.log.")
log_info("rMATS done for all comparisons. Next: Rscript splicing/scripts/run_splicing.R --config <comparison>/splicing_config.yaml")
