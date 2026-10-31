#!/usr/bin/env Rscript
# =============================================================================
# run_isoforms.R  --  differential isoform expression and usage (one run)
# =============================================================================
# The entry point of the isoform module, built like scripts/run_pipeline.R:
#   1. creates a run folder  runs/<YYYYmmdd-HHMMSS>_<run_name>/
#   2. merges parameters (isoform defaults < --config < --based_on run < flags)
#   3. runs I01-I04, each in its own R process with its own log (I01 copies the
#      quantification files and the sample table into the run)
#   4. writes provenance.json and run_summary.json
#
# Usage
#   Rscript isoforms/scripts/run_isoforms.R --quant_dir results/salmon --metadata_file samples.csv \
#       --gtf_file ref/genes.gtf --condition genotype --reference_level WT --test_level KO
#   Rscript isoforms/scripts/run_isoforms.R --based_on runs/<run> --set "min_delta_usage=0.2"
#
# Input: Salmon (or kallisto) output, one folder per sample (later/scripts/U02_salmon_quant.R makes it).
# Exit code: 0 = critical steps passed, 1 = a critical step failed.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "isoforms/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "R", "gtf.R"))
source(file.path(REPO_ROOT, "isoforms", "R", "isoform_utils.R"))
suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config", type = "character", help = "Run config YAML (values that differ from isoforms/config/default_params.yaml)"),
  make_option("--based_on", type = "character", help = "Start from the parameters of an earlier isoform run folder"),
  make_option("--runs_dir", type = "character", help = "Where run folders are created (default: <repo>/runs)"),
  make_option("--run_name", type = "character", help = "Short label for the run folder name"),
  make_option("--steps", type = "character", default = "1-4", help = "Steps to run, e.g. 1-4 [default %default]"),
  make_option("--set", type = "character", help = "Overrides, e.g. \"min_delta_usage=0.2;covariates=[batch]\""),
  make_option("--quant_dir", type = "character", help = "Folder with <sample>/quant.sf"),
  make_option("--metadata_file", type = "character", help = "Sample table (first column = sample IDs)"),
  make_option("--gtf_file", type = "character", help = "Annotation GTF"),
  make_option("--condition", type = "character", help = "Sample-table column with the groups"),
  make_option("--reference_level", type = "character", help = "Control / baseline group"),
  make_option("--test_level", type = "character", help = "Group compared to the reference")
), description = "Differential isoform expression and usage."))

steps <- local({
  s <- gsub(" ", "", opt$steps)
  n <- unlist(lapply(strsplit(s, ",")[[1]], function(p) if (grepl("-", p)) { r <- as.integer(strsplit(p, "-")[[1]]); seq(r[1], r[2]) } else as.integer(p)))
  n <- sort(unique(n))
  if (any(is.na(n)) || any(n < 1 | n > length(ISOFORM_STEPS))) stop("--steps must use numbers 1 to ", length(ISOFORM_STEPS))
  n
})

# ---- Parameters -------------------------------------------------------------------
defaults <- load_isoform_defaults()
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
cli <- opt[intersect(c("quant_dir", "metadata_file", "gtf_file", "condition", "reference_level", "test_level", "run_name"), names(opt))]
params <- modifyList(params, c(cli[!vapply(cli, is.null, logical(1))], parse_set_string(opt$set)))
for (k in ISOFORM_LIST_PARAMS) if (k %in% names(params)) params[[k]] <- as_list_param(params[[k]])
unknown <- setdiff(names(params), c(names(defaults), "based_on"))
if (length(unknown)) log_warn("Unknown parameter(s) ignored (typo?): ", paste(unknown, collapse = ", "))

# ---- Run folder ----------------------------------------------------------------------------
for (k in c("quant_dir", "metadata_file", "gtf_file", "tx2gene_file")) {
  p <- resolve_path(params[[k]])
  if (!is.null(p) && file.exists(p)) params[[k]] <- normalizePath(p)
}
runs_dir <- opt$runs_dir %||% file.path(REPO_ROOT, "runs")
dir.create(runs_dir, showWarnings = FALSE, recursive = TRUE)
run_dir <- file.path(normalizePath(runs_dir), paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "_", make.names(params$run_name %||% "isoforms")))
while (dir.exists(run_dir)) { Sys.sleep(1); run_dir <- paste0(run_dir, "b") }
dir.create(file.path(run_dir, "inputs"), recursive = TRUE)
dir.create(file.path(run_dir, "logs"))
if (!is.null(params$metadata_file) && file.exists(params$metadata_file)) {
  file.copy(params$metadata_file, file.path(run_dir, "inputs", basename(params$metadata_file)))
}
write_yaml_file(params, file.path(run_dir, "params.yaml"),
                header = c("Parameters requested for this isoform run (defaults + config + overrides).",
                           "Step I01 fills in the rest: see I01_validate_isoforms/resolved_params.yaml"))
log_info("New isoform run: ", run_dir)

# ---- Steps --------------------------------------------------------------------------
rscript <- file.path(R.home("bin"), "Rscript")
records <- list()
overall <- "success"
for (i in steps) {
  st <- ISOFORM_STEPS[[i]]
  log_file <- file.path(run_dir, "logs", paste0(st$name, ".log"))
  log_info(sprintf("Step %d/%d: %s ...", i, length(ISOFORM_STEPS), st$name))
  t0 <- Sys.time()
  code <- system(paste(shQuote(rscript), shQuote(file.path(ISOFORM_ROOT, "scripts", paste0(st$name, ".R"))), "--run_dir", shQuote(run_dir),
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
code_bank <- list(version = CODE_BANK_VERSION, module = "isoforms",
                  commit = if (length(commit) && is.null(attr(commit, "status"))) commit[1] else NULL,
                  uncommitted_changes = length(git("status", "--porcelain", "--untracked-files=no")) > 0)
resolved_path <- file.path(run_dir, "I01_validate_isoforms", "resolved_params.yaml")
all_pkgs <- list()
for (s in names(records)) {
  p <- read_json_file(file.path(run_dir, "provenance", paste0(s, ".json")))
  if (!is.null(p$packages)) all_pkgs <- modifyList(all_pkgs, as.list(p$packages))
}
write_json_file(list(
  run_id = basename(run_dir), module = "isoforms", status = overall, finished = timestamp_now(),
  code_bank = code_bank, r_version = R.version.string, platform = R.version$platform, packages = all_pkgs,
  parameters_requested = yaml::read_yaml(file.path(run_dir, "params.yaml")),
  parameters_resolved = if (file.exists(resolved_path)) yaml::read_yaml(resolved_path) else NULL,
  input_md5 = md5_of(list.files(file.path(run_dir, "inputs"), full.names = TRUE, recursive = TRUE)),
  steps = records,
  rerun_command = sprintf("Rscript isoforms/scripts/run_isoforms.R --based_on runs/%s", basename(run_dir))
), file.path(run_dir, "provenance.json"))

val <- read_json_file(file.path(run_dir, "I01_validate_isoforms", "validation_report.json"))
dei <- read_json_file(file.path(run_dir, "I02_differential_isoform_expression", "dei_summary.json"))
dui <- read_json_file(file.path(run_dir, "I03_differential_isoform_usage", "dui_summary.json"))
warn <- if (!is.null(val) && length(val$checks)) {
  ch <- val$checks
  ch <- ch[ch$status %in% c("warning", "fail"), , drop = FALSE]
  if (NROW(ch)) as.list(paste0("[validation:", ch$id, "] ", ch$message)) else list()
} else list()
write_json_file(list(
  run_id = basename(run_dir), module = "isoforms", status = overall, steps = lapply(records, `[[`, "status"),
  comparison = dei$comparison, design_formula = dei$design_formula,
  dei = dei[intersect(c("significance_rule", "n_transcripts_tested", "n_significant", "n_up", "n_down", "n_genes_with_significant_transcripts"), names(dei))],
  dui = dui[intersect(c("status", "reason", "significance_rule", "n_genes_tested", "n_genes_significant", "n_isoform_switches", "n_transcripts_significant"), names(dui))],
  warnings = warn
), file.path(run_dir, "run_summary.json"))
log_info("Isoform run finished: ", overall, " -> ", run_dir)
quit(save = "no", status = if (overall == "failed") 1 else 0)
