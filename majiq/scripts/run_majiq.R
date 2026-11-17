#!/usr/bin/env Rscript
# =============================================================================
# run_majiq.R  --  differential splicing with MAJIQ v3 (one run)
# =============================================================================
# The entry point of the MAJIQ module, built like scripts/run_pipeline.R:
#   1. creates a run folder  runs/<YYYYmmdd-HHMMSS>_<run_name>/
#   2. merges parameters (MAJIQ defaults < --config < --based_on run < flags)
#   3. runs MJ01-MJ06, each in its own R process with its own log
#        MJ01 check the inputs, MAJIQ and its licence; BAMs linked as <sample>.bam
#        MJ02 majiq build (splice graph from the GFF3 and the BAMs)
#        MJ03 majiq psi-coverage per group (+ majiq-v3 sg-coverage for modulize)
#        MJ04 majiq deltapsi; voila tsv (all and changing); voila modulize
#        MJ05 tables, summary and figures
#        MJ06 methods paragraph
#   4. writes provenance.json and run_summary.json
#
# Usage
#   Rscript majiq/scripts/run_majiq.R --bam_dir results/star --metadata_file samples.csv \
#       --gff3_file ref/genes.gff3.gz --condition genotype --reference_level WT --test_level KO
#   Rscript majiq/scripts/run_majiq.R --based_on runs/<run> --set "probability_threshold=0.9"
#
# MAJIQ itself (academic licence): majiq/install/install_majiq.sh.
# Exit code: 0 = critical steps passed, 1 = a critical step failed.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "majiq/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "majiq", "R", "majiq_utils.R"))
suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config", type = "character", help = "Run config YAML (values that differ from majiq/config/default_params.yaml)"),
  make_option("--based_on", type = "character", help = "Start from the parameters of an earlier MAJIQ run folder"),
  make_option("--runs_dir", type = "character", help = "Where run folders are created (default: <repo>/runs)"),
  make_option("--run_name", type = "character", help = "Short label for the run folder name"),
  make_option("--steps", type = "character", default = "1-6", help = "Steps to run, e.g. 1-6 [default %default]"),
  make_option("--set", type = "character", help = "Overrides, e.g. \"probability_threshold=0.9;threads=16\""),
  make_option("--bam_dir", type = "character", help = "Folder with the sorted, indexed BAMs"),
  make_option("--metadata_file", type = "character", help = "Sample table (first column = sample IDs)"),
  make_option("--gff3_file", type = "character", help = "GFF3 annotation (same release as the GTF used for alignment)"),
  make_option("--condition", type = "character", help = "Sample-table column with the groups"),
  make_option("--reference_level", type = "character", help = "Control / baseline group"),
  make_option("--test_level", type = "character", help = "Group compared to the reference")
), description = "Differential splicing with MAJIQ v3."))

steps <- local({
  s <- gsub(" ", "", opt$steps)
  n <- unlist(lapply(strsplit(s, ",")[[1]], function(p) if (grepl("-", p)) { r <- as.integer(strsplit(p, "-")[[1]]); seq(r[1], r[2]) } else as.integer(p)))
  n <- sort(unique(n))
  if (any(is.na(n)) || any(n < 1 | n > length(MAJIQ_STEPS))) stop("--steps must use numbers 1 to ", length(MAJIQ_STEPS))
  n
})

# ---- Parameters -------------------------------------------------------------------
defaults <- load_majiq_defaults()
params <- defaults
if (!is.null(opt$based_on)) {
  base <- file.path(opt$based_on, "params.yaml")
  if (!file.exists(base)) stop("--based_on folder has no params.yaml: ", opt$based_on)
  params <- modifyList(params, yaml::read_yaml(base))
  params$based_on <- basename(normalizePath(opt$based_on))
}
if (!is.null(opt$config)) {
  cfg <- resolve_path(opt$config)
  if (!file.exists(cfg)) stop("Config file not found: ", opt$config)
  params <- modifyList(params, yaml::read_yaml(cfg) %||% list())
}
cli <- opt[intersect(c("bam_dir", "metadata_file", "gff3_file", "condition", "reference_level", "test_level", "run_name"), names(opt))]
params <- modifyList(params, c(cli[!vapply(cli, is.null, logical(1))], parse_set_string(opt$set)))
for (k in MAJIQ_LIST_PARAMS) if (k %in% names(params)) params[[k]] <- as_list_param(params[[k]])
unknown <- setdiff(names(params), c(names(defaults), "based_on"))
if (length(unknown)) log_warn("Unknown parameter(s) ignored (typo?): ", paste(unknown, collapse = ", "))

# ---- Run folder ----------------------------------------------------------------------------
for (k in c("bam_dir", "metadata_file", "gff3_file")) {
  p <- resolve_path(params[[k]])
  if (!is.null(p) && file.exists(p)) params[[k]] <- normalizePath(p)
}
runs_dir <- opt$runs_dir %||% file.path(REPO_ROOT, "runs")
dir.create(runs_dir, showWarnings = FALSE, recursive = TRUE)
run_dir <- file.path(normalizePath(runs_dir), paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "_", make.names(params$run_name %||% "majiq")))
while (dir.exists(run_dir)) { Sys.sleep(1); run_dir <- paste0(run_dir, "b") }
dir.create(file.path(run_dir, "inputs"), recursive = TRUE)
dir.create(file.path(run_dir, "logs"))
if (!is.null(params$metadata_file) && file.exists(params$metadata_file)) {
  file.copy(params$metadata_file, file.path(run_dir, "inputs", basename(params$metadata_file)))
}
write_yaml_file(params, file.path(run_dir, "params.yaml"),
                header = c("Parameters requested for this MAJIQ run (defaults + config + overrides).",
                           "Step MJ01 fills in the rest: see MJ01_validate_majiq/resolved_params.yaml"))
log_info("New MAJIQ run: ", run_dir)

# ---- Steps --------------------------------------------------------------------------
rscript <- file.path(R.home("bin"), "Rscript")
records <- list()
overall <- "success"
for (i in steps) {
  st <- MAJIQ_STEPS[[i]]
  log_file <- file.path(run_dir, "logs", paste0(st$name, ".log"))
  log_info(sprintf("Step %d/%d: %s ...", i, length(MAJIQ_STEPS), st$name))
  t0 <- Sys.time()
  code <- system(paste(shQuote(rscript), shQuote(file.path(MAJIQ_ROOT, "scripts", paste0(st$name, ".R"))), "--run_dir", shQuote(run_dir),
                       ">", shQuote(log_file), "2>&1"))
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  ok <- identical(as.integer(code), 0L)
  records[[st$name]] <- list(status = if (ok) "success" else "failed", exit_code = code, seconds = secs,
                             critical = st$critical, log = file.path("logs", basename(log_file)))
  if (ok) {
    log_info(sprintf("  done in %.1f s", secs))
  } else {
    log_msg(sprintf("  FAILED (exit code %s). Last lines of %s:", code, log_file), level = "ERROR")
    cat(paste0("    | ", tail(readLines(log_file, warn = FALSE), 15), "\n"), sep = "")
    if (st$critical) { overall <- "failed"; break }
    overall <- "success_with_errors"
  }
}

# ---- Provenance and summary ------------------------------------------------------------
git <- function(...) tryCatch(suppressWarnings(system2("git", c("-C", shQuote(REPO_ROOT), ...), stdout = TRUE, stderr = FALSE)),
                              error = function(e) character(0))
commit <- git("rev-parse", "HEAD")
code_bank <- list(version = CODE_BANK_VERSION, module = "majiq",
                  commit = if (length(commit) && is.null(attr(commit, "status"))) commit[1] else NULL,
                  uncommitted_changes = length(git("status", "--porcelain", "--untracked-files=no")) > 0)
resolved_path <- file.path(run_dir, "MJ01_validate_majiq", "resolved_params.yaml")
all_pkgs <- list()
for (s in names(records)) {
  p <- read_json_file(file.path(run_dir, "provenance", paste0(s, ".json")))
  if (!is.null(p$packages)) all_pkgs <- modifyList(all_pkgs, as.list(p$packages))
}
write_json_file(list(
  run_id = basename(run_dir), module = "majiq", status = overall, finished = timestamp_now(),
  code_bank = code_bank, r_version = R.version.string, platform = R.version$platform, packages = all_pkgs,
  parameters_requested = yaml::read_yaml(file.path(run_dir, "params.yaml")),
  parameters_resolved = if (file.exists(resolved_path)) yaml::read_yaml(resolved_path) else NULL,
  input_md5 = md5_of(list.files(file.path(run_dir, "inputs"), full.names = TRUE, recursive = TRUE)),
  steps = records,
  rerun_command = sprintf("Rscript majiq/scripts/run_majiq.R --based_on runs/%s", basename(run_dir))
), file.path(run_dir, "provenance.json"))

val <- read_json_file(file.path(run_dir, "MJ01_validate_majiq", "validation_report.json"))
sm <- read_json_file(file.path(run_dir, "MJ05_majiq_results", "majiq_summary.json"))
warn <- if (!is.null(val) && length(val$checks)) {
  ch <- val$checks
  ch <- ch[ch$status %in% c("warning", "fail"), , drop = FALSE]
  if (NROW(ch)) as.list(paste0("[validation:", ch$id, "] ", ch$message)) else list()
} else list()
write_json_file(list(
  run_id = basename(run_dir), module = "majiq", status = overall, steps = lapply(records, `[[`, "status"),
  comparison = sm$comparison,
  majiq = sm[intersect(c("significance_rule", "sign_convention", "n_lsvs_quantified", "n_lsvs_changing", "n_genes_changing",
                         "modulize_by_type"), names(sm))],
  warnings = warn
), file.path(run_dir, "run_summary.json"))
log_info("MAJIQ run finished: ", overall, " -> ", run_dir)
quit(save = "no", status = if (overall == "failed") 1 else 0)
