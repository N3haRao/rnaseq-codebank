#!/usr/bin/env Rscript
# =============================================================================
# SU01_validate_suppa.R  --  check the inputs and make SUPPA's TPM tables
# =============================================================================
# Checks the sample table and the two groups, finds each sample's Salmon quant.sf,
# checks that the transcripts match the GTF, finds SUPPA2, and writes one TPM table
# per group in SUPPA's format (every GTF transcript; transcripts Salmon did not
# report get 0, as in the lab's scripts).
#
# Outputs (in <run_dir>/SU01_validate_suppa/): validation_report.json / .md,
# samples_validated.csv, <reference>.tpm, <test>.tpm, transcripts.csv
# (transcript -> gene, gene name), annotation.gtf (only when the GTF was .gz),
# resolved_params.yaml
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "suppa/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "R", "gtf.R"))
source(file.path(REPO_ROOT, "suppa", "R", "suppa_utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "SU01_validate_suppa"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "SUPPA step SU01: check the inputs."))
if (is_blank(opt$run_dir)) fail("--run_dir is required.")
params <- load_suppa_params(opt$run_dir, use_resolved = FALSE)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
log_info("Step SU01: checking the inputs for SUPPA2")

checks <- list()
check <- function(id, status, message) {
  checks[[length(checks) + 1]] <<- list(id = id, status = status, message = message)
  (if (status == "pass") log_info else log_warn)(sprintf("[%s] %s: %s", status, id, message))
}
finish <- function() {
  st <- vapply(checks, `[[`, "", "status")
  overall <- if (any(st == "fail")) "fail" else if (any(st == "warning")) "warning" else "pass"
  write_json_file(list(status = overall, n_fail = sum(st == "fail"), n_warning = sum(st == "warning"), checks = checks),
                  file.path(out_dir, "validation_report.json"))
  writeLines(c("# SUPPA2 input check", "", sprintf("Overall: **%s**", overall), "",
               vapply(checks, function(c) sprintf("- **%s** (%s): %s", c$id, c$status, c$message), "")),
             file.path(out_dir, "validation_report.md"))
  overall
}
stop_now <- function(msg) { finish(); fail(msg) }

# ---- SUPPA2 itself ----------------------------------------------------------------------
cmd <- find_suppa(params)
if (is.null(cmd)) {
  check("suppa", "fail", paste("SUPPA2 was not found. Install it with suppa/install/install_suppa.sh, or set suppa_command",
                               "to the suppa.py to use."))
  stop_now("SUPPA2 not found.")
}
ver_log <- file.path(out_dir, "suppa_version.log")
code <- suppressWarnings(system2(cmd[1], c(cmd[-1], "-h"), stdout = ver_log, stderr = ver_log))
if (!identical(as.integer(code), 0L)) {
  check("suppa", "fail", sprintf("SUPPA2 (%s) did not start; see %s.", paste(cmd, collapse = " "), basename(ver_log)))
  stop_now("SUPPA2 does not run.")
}
check("suppa", "pass", sprintf("SUPPA2: %s", paste(cmd, collapse = " ")))

# ---- Sample table and groups ------------------------------------------------------------
mf <- resolve_path(params$metadata_file)
if (is.null(mf) || !file.exists(mf)) { check("metadata_file", "fail", "The sample table was not found."); stop_now("No sample table.") }
meta <- read_table_file(mf)
names(meta)[1] <- "sample_id"
meta[] <- lapply(meta, function(x) trimws(as.character(x)))
excl <- intersect(params$exclude_samples, meta$sample_id)
if (length(excl)) { meta <- meta[!meta$sample_id %in% excl, ]; check("excluded_samples", "pass", paste("Left out:", paste(excl, collapse = ", "))) }
if (anyDuplicated(meta$sample_id)) { check("duplicate_samples", "fail", "Sample IDs are repeated in the sample table."); stop_now("Duplicate samples.") }
cond <- params$condition
if (!cond %in% names(meta)) { check("condition_column", "fail", sprintf("No column '%s' in the sample table (columns: %s).", cond, paste(names(meta)[-1], collapse = ", "))); stop_now("No condition column.") }
lv <- unique(meta[[cond]])
ref <- params$reference_level; test <- params$test_level
if (is_blank(ref) || !ref %in% lv) { check("reference_level", "fail", sprintf("reference_level must be one of: %s.", paste(lv, collapse = ", "))); stop_now("Bad reference level.") }
if (is_blank(test)) test <- if (length(lv) == 2) setdiff(lv, ref) else NULL
if (is.null(test) || !test %in% lv || test == ref) { check("test_level", "fail", sprintf("test_level must be another value of '%s': %s.", cond, paste(setdiff(lv, ref), collapse = ", "))); stop_now("Bad test level.") }
meta <- meta[meta[[cond]] %in% c(ref, test), , drop = FALSE]
n_ref <- sum(meta[[cond]] == ref); n_test <- sum(meta[[cond]] == test)
if (min(n_ref, n_test) < 2) { check("replicates", "fail", sprintf("%s has %d and %s has %d sample(s); SUPPA2's test needs at least 2 per group.", ref, n_ref, test, n_test)); stop_now("Too few replicates.") }
check("replicates", if (min(n_ref, n_test) < 3) "warning" else "pass",
      sprintf("%s: %d samples, %s: %d samples%s.", ref, n_ref, test, n_test, if (min(n_ref, n_test) < 3) "; at least 3 per group are recommended" else ""))
bad_types <- setdiff(params$event_types, c("SE", "SS", "MX", "RI", "FL"))
if (!length(params$event_types) || length(bad_types)) { check("event_types", "fail", "event_types must be some of SE, SS, MX, RI, FL."); stop_now("Bad event types.") }

# ---- Annotation -------------------------------------------------------------------------
gtf <- resolve_path(params$gtf_file)
if (is.null(gtf) || !file.exists(gtf)) { check("gtf_file", "fail", "gtf_file was not found: SUPPA2 finds the events in the annotation."); stop_now("No GTF.") }
if (grepl("\\.gz$", gtf)) {   # SUPPA reads plain text
  plain <- file.path(out_dir, "annotation.gtf")
  zin <- gzfile(gtf, "r"); zout <- file(plain, "w")
  repeat { x <- readLines(zin, n = 200000, warn = FALSE); if (!length(x)) break; writeLines(x, zout) }
  close(zin); close(zout)
  gtf <- plain
}
tx <- tx2gene_from_gtf(gtf)
if (!nrow(tx)) { check("gtf_file", "fail", "No transcripts (transcript_id) in the GTF."); stop_now("Empty GTF.") }
check("gtf_file", "pass", sprintf("%d transcripts in %d genes.", nrow(tx), length(unique(tx$gene_id))))
write_table(tx, file.path(out_dir, "transcripts.csv"))

# ---- Salmon quantification ----------------------------------------------------------------
qdir <- resolve_path(params$quant_dir)
if (is.null(qdir) || !dir.exists(qdir)) { check("quant_dir", "fail", "quant_dir was not found."); stop_now("No quantification folder.") }
files <- setNames(file.path(qdir, meta$sample_id, "quant.sf"), meta$sample_id)
if (!all(file.exists(files))) {
  check("quant_files", "fail", sprintf("No quant.sf for sample(s): %s (expected <quant_dir>/<sample>/quant.sf).",
                                       paste(names(files)[!file.exists(files)], collapse = ", ")))
  stop_now("Missing quantification files.")
}
check("quant_files", "pass", sprintf("quant.sf found for all %d samples.", length(files)))

read_tpm <- function(f) {
  q <- data.table::fread(f, data.table = FALSE, select = c("Name", "TPM"))
  setNames(q$TPM, sub("\\|.*$", "", q$Name))            # GENCODE transcript FASTA names: ENST...|ENSG...|...
}
tpm <- lapply(files, read_tpm)
ids <- names(tpm[[1]])
# transcripts named as in the GTF? (with or without the .version)
matched <- mean(ids %in% tx$transcript_id)
strip_ver <- function(x) sub("\\.\\d+$", "", x)
map_ids <- NULL
if (matched < 0.9) {
  m2 <- mean(strip_ver(ids) %in% strip_ver(tx$transcript_id))
  if (m2 >= 0.9) {
    map_ids <- setNames(tx$transcript_id, strip_ver(tx$transcript_id))
    check("transcript_ids", "warning", sprintf("Transcript versions differ between Salmon and the GTF; matched without versions (%.1f%%).", 100 * m2))
    matched <- m2
  } else {
    check("transcript_ids", "fail", sprintf("Only %.0f%% of the quantified transcripts are in the GTF. Use the GTF the Salmon index was built from.", 100 * max(matched, m2)))
    stop_now("Annotation does not match the quantification.")
  }
} else check("transcript_ids", "pass", sprintf("%.1f%% of %d quantified transcripts are in the GTF.", 100 * matched, length(ids)))

# one matrix: GTF transcripts x samples, 0 where Salmon has no value
M <- matrix(0, nrow = nrow(tx), ncol = length(files), dimnames = list(tx$transcript_id, names(files)))
for (s in names(tpm)) {
  v <- tpm[[s]]
  nm <- if (is.null(map_ids)) names(v) else unname(map_ids[strip_ver(names(v))])
  keep <- !is.na(nm) & nm %in% rownames(M)
  M[nm[keep], s] <- v[keep]
}
lib_low <- names(which(colSums(M > 1) < 1000))
if (length(lib_low)) check("expressed_transcripts", "warning", sprintf("Fewer than 1,000 transcripts above 1 TPM in: %s.", paste(lib_low, collapse = ", ")))
for (g in c(ref, test)) {
  s <- meta$sample_id[meta[[cond]] == g]
  write_suppa_matrix(M[, s, drop = FALSE], file.path(out_dir, paste0(make.names(g), ".tpm")))
}
write_table(meta, file.path(out_dir, "samples_validated.csv"))

# ---- Optional rMATS tables ------------------------------------------------------------------
rd <- resolve_path(params$rmats_dir)
if (!is_blank(params$rmats_dir)) {
  ok <- !is.null(rd) && dir.exists(rd) && length(list.files(rd, pattern = sprintf("\\.MATS\\.%s\\.txt$", params$rmats_count_type)))
  check("rmats_dir", if (ok) "pass" else "warning",
        if (ok) "rMATS tables found: step SU04 compares SUPPA2 with rMATS." else
          sprintf("No *.MATS.%s.txt in rmats_dir: the comparison with rMATS (SU04) is skipped.", params$rmats_count_type))
}

organism <- params$organism %||% "auto"
if (organism == "auto") organism <- if (mean(grepl("^ENSMUSG", head(tx$gene_id, 500))) > 0.5) "mouse" else
  if (mean(grepl("^ENSG\\d", head(tx$gene_id, 500))) > 0.5) "human" else "other"

if (finish() == "fail") fail("The input check failed.")
write_yaml_file(list(reference_level = ref, test_level = test, organism = organism, gtf_used = gtf,
                     suppa_command = paste(cmd, collapse = " "),
                     tpm_files = list(reference = paste0(make.names(ref), ".tpm"), test = paste0(make.names(test), ".tpm")),
                     n_samples = nrow(meta), n_transcripts = nrow(tx),
                     group_sizes = setNames(list(n_ref, n_test), c(ref, test))),
                file.path(out_dir, "resolved_params.yaml"), header = "Values used by all later SUPPA steps (filled in by SU01).")
write_step_provenance(run_dir, STEP, started, params_used = params[c("condition", "event_types", "exclude_samples")],
                      inputs = unname(files), outputs = list.files(out_dir))
log_info(sprintf("Step SU01 done: %d samples (%s %d, %s %d), %d transcripts.", nrow(meta), ref, n_ref, test, n_test, nrow(tx)))
