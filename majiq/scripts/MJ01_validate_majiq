#!/usr/bin/env Rscript
# =============================================================================
# MJ01_validate_majiq.R  --  check MAJIQ, its licence and the inputs
# =============================================================================
# Finds majiq / voila and the licence file, checks the sample table and the two
# groups, finds each sample's sorted BAM and its index, and links them into the
# run as bams/<sample>.bam (MAJIQ names its files after the BAMs), indexing with
# samtools when there is no .bai. Writes MAJIQ's experiments.tsv (path, group).
#
# Outputs (in <run_dir>/MJ01_validate_majiq/): validation_report.json / .md,
# samples_validated.csv, bams/, experiments.tsv, majiq_version.log, resolved_params.yaml
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "majiq/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "majiq", "R", "majiq_utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "MJ01_validate_majiq"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "MAJIQ step MJ01: check the inputs."))
if (is_blank(opt$run_dir)) fail("--run_dir is required.")
params <- load_majiq_params(opt$run_dir, use_resolved = FALSE)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
log_info("Step MJ01: checking MAJIQ and the inputs")

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
  writeLines(c("# MAJIQ input check", "", sprintf("Overall: **%s**", overall), "",
               vapply(checks, function(c) sprintf("- **%s** (%s): %s", c$id, c$status, c$message), "")),
             file.path(out_dir, "validation_report.md"))
  overall
}
stop_now <- function(msg) { finish(); fail(msg) }

# ---- MAJIQ, voila and the licence -----------------------------------------------------------
bin <- find_majiq_bin(params)
if (is.null(bin)) {
  check("majiq", "fail", paste("MAJIQ was not found. Install it with majiq/install/install_majiq.sh (academic licence),",
                               "or set majiq_bin_dir to the folder holding majiq and voila."))
  stop_now("MAJIQ not found.")
}
need <- c("majiq", "voila", if (isTRUE(params$run_modulize)) "majiq-v3")
miss <- need[!file.exists(file.path(bin, need))]
if (length(miss)) { check("majiq", "fail", sprintf("Not in %s: %s.", bin, paste(miss, collapse = ", "))); stop_now("MAJIQ incomplete.") }
lic <- find_majiq_license(params)
if (is.null(lic)) {
  check("majiq_license", "fail", paste("No MAJIQ licence file. MAJIQ >= 2.5 does not run without one: put majiq_license*.lic in",
                                       "$COPILOT_HOME (see majiq/install/install_majiq.sh) or set majiq_license_file."))
  stop_now("No MAJIQ licence.")
}
check("majiq_license", "pass", paste("Licence file:", basename(lic)))
ver_log <- file.path(out_dir, "majiq_version.log")
code <- suppressWarnings(system2(file.path(bin, "majiq"), "--version", stdout = ver_log, stderr = ver_log))
ver <- if (file.exists(ver_log)) grep("[0-9]+\\.[0-9]+", readLines(ver_log, warn = FALSE), value = TRUE) else character(0)
ver <- if (length(ver)) regmatches(ver[length(ver)], regexpr("[0-9]+\\.[0-9]+(\\.[0-9]+)?", ver[length(ver)])) else "unknown"
if (!identical(as.integer(code), 0L)) { check("majiq", "fail", sprintf("majiq --version failed; see %s.", basename(ver_log))); stop_now("MAJIQ does not run.") }
check("majiq", "pass", sprintf("MAJIQ %s in %s.", ver, bin))

# ---- Sample table and groups ------------------------------------------------------------------
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
if (min(n_ref, n_test) < 2) { check("replicates", "fail", sprintf("%s has %d and %s has %d sample(s); both need at least 2.", ref, n_ref, test, n_test)); stop_now("Too few replicates.") }
check("replicates", if (min(n_ref, n_test) < 3) "warning" else "pass",
      sprintf("%s: %d samples, %s: %d samples%s.", ref, n_ref, test, n_test, if (min(n_ref, n_test) < 3) "; at least 3 per group are recommended" else ""))
if (!params$strandness %in% c("AUTO", "FORWARD", "REVERSE", "NONE")) { check("strandness", "fail", "strandness must be AUTO, FORWARD, REVERSE or NONE."); stop_now("Bad strandness.") }

# ---- Annotation ---------------------------------------------------------------------------------
gff <- resolve_path(params$gff3_file)
if (is.null(gff) || !file.exists(gff)) { check("gff3_file", "fail", "gff3_file was not found: MAJIQ needs the annotation as GFF3 (e.g. gencode.vXX.primary_assembly.annotation.gff3.gz)."); stop_now("No GFF3.") }
first <- readLines(if (grepl("\\.gz$", gff)) gzfile(gff) else gff, n = 200, warn = FALSE)
if (!any(grepl("\tgene\t|\tmRNA\t|\ttranscript\t", first)) || !any(grepl("ID=|Parent=", first))) {
  check("gff3_file", "fail", "The annotation does not look like GFF3 (no ID= / Parent= attributes). A GTF will not work; use the GFF3 of the same release.")
  stop_now("Not a GFF3.")
}
check("gff3_file", "pass", basename(gff))

# ---- BAMs ---------------------------------------------------------------------------------------
bdir <- resolve_path(params$bam_dir)
if (is.null(bdir) || !dir.exists(bdir)) { check("bam_dir", "fail", "bam_dir was not found."); stop_now("No BAM folder.") }
bams <- vapply(meta$sample_id, function(s) find_sample_bam(bdir, s), "")
if (any(is.na(bams))) {
  check("bam_files", "fail", sprintf("No BAM for: %s (expected <sample>.bam, <sample>.sorted.bam or <sample>/Aligned.sortedByCoord.out.bam).",
                                     paste(meta$sample_id[is.na(bams)], collapse = ", ")))
  stop_now("Missing BAMs.")
}
link_dir <- file.path(out_dir, "bams")
dir.create(link_dir, showWarnings = FALSE)
samtools <- Sys.which("samtools")
no_index <- character(0)
for (i in seq_along(bams)) {
  s <- meta$sample_id[i]
  dest <- file.path(link_dir, paste0(s, ".bam"))
  unlink(c(dest, paste0(dest, ".bai")))
  file.symlink(normalizePath(bams[i]), dest)
  bai <- c(paste0(bams[i], ".bai"), sub("\\.bam$", ".bai", bams[i]))
  bai <- bai[file.exists(bai)]
  if (length(bai)) file.symlink(normalizePath(bai[1]), paste0(dest, ".bai"))
  else if (nzchar(samtools)) {
    log_info("Indexing ", s, " (samtools index)")
    code <- system2(samtools, c("index", "-@", "4", "-o", paste0(dest, ".bai"), dest))
    if (code != 0) no_index <- c(no_index, s)
  } else no_index <- c(no_index, s)
}
if (length(no_index)) { check("bam_index", "fail", sprintf("No BAM index (.bai) for: %s, and samtools could not make one. Run samtools index on the BAMs.", paste(no_index, collapse = ", "))); stop_now("BAMs not indexed.") }
check("bam_files", "pass", sprintf("Sorted, indexed BAMs found for all %d samples.", nrow(meta)))
grp <- make.names(meta[[cond]])
write.table(data.frame(path = file.path(link_dir, paste0(meta$sample_id, ".bam")), group = grp),
            file.path(out_dir, "experiments.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write_table(meta, file.path(out_dir, "samples_validated.csv"))

if (finish() == "fail") fail("The input check failed.")
organism <- params$organism %||% "auto"
if (organism == "auto") {
  g <- grep("gene_id=|ID=", first, value = TRUE)
  organism <- if (any(grepl("ENSMUSG", g))) "mouse" else if (any(grepl("ENSG\\d", g))) "human" else "other"
}
write_yaml_file(list(reference_level = ref, test_level = test, ref_group = make.names(ref), test_group = make.names(test),
                     organism = organism, majiq_bin_dir = bin, majiq_license_file = lic, majiq_version = ver, gff3_used = gff,
                     samples = list(reference = meta$sample_id[meta[[cond]] == ref], test = meta$sample_id[meta[[cond]] == test]),
                     n_samples = nrow(meta), group_sizes = setNames(list(n_ref, n_test), c(ref, test))),
                file.path(out_dir, "resolved_params.yaml"), header = "Values used by all later MAJIQ steps (filled in by MJ01).")
write_step_provenance(run_dir, STEP, started, params_used = params[c("condition", "strandness", "exclude_samples")],
                      outputs = list.files(out_dir), notes = list(majiq_version = ver))
log_info(sprintf("Step MJ01 done: %d samples (%s %d, %s %d), MAJIQ %s.", nrow(meta), ref, n_ref, test, n_test, ver))
