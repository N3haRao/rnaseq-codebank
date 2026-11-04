#!/usr/bin/env Rscript
# =============================================================================
# run_compare.R  --  compare the results of two or more RNA-seq studies
# =============================================================================
# The entry point of the comparison module, built like scripts/run_pipeline.R:
#   1. creates a run folder  runs/<YYYYmmdd-HHMMSS>_<run_name>/
#   2. merges parameters (compare defaults < --config < --based_on run < flags)
#      and saves them as params.yaml
#   3. copies each study's results into inputs/<study>/ so the run is self-contained
#   4. runs C01-C05, each in its own R process with its own log
#   5. writes provenance.json and run_summary.json
#
# Usage
#   Rscript compare/scripts/run_compare.R --studies "ours=runs/<run1>,GSE60450=runs/<run2>"
#   Rscript compare/scripts/run_compare.R --config my_comparison.yaml
#   Rscript compare/scripts/run_compare.R --based_on runs/<compare_run> --set "padj_cutoff=0.1"
#
# Studies are compared by their results (fold changes, significant genes, gene
# rankings, pathways), never by merging their raw counts.
# Exit code: 0 = critical steps passed, 1 = a critical step failed.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "compare/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "compare", "R", "compare_utils.R"))
suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config", type = "character", help = "Run config YAML (values that differ from compare/config/default_params.yaml)"),
  make_option("--based_on", type = "character", help = "Start from the parameters of an earlier comparison run folder"),
  make_option("--runs_dir", type = "character", help = "Where run folders are created (default: <repo>/runs)"),
  make_option("--run_name", type = "character", help = "Short label for the run folder name"),
  make_option("--steps", type = "character", default = "1-5", help = "Steps to run, e.g. 1-5 or 2-5 [default %default]"),
  make_option("--set", type = "character", help = "Overrides, e.g. \"padj_cutoff=0.1;rank_by=signed_p\""),
  make_option("--studies", type = "character",
              help = "label=path pairs separated by commas; path = a pipeline run folder or a DE results table")
), description = "Compare differential-expression results across studies."))

steps <- local({
  s <- gsub(" ", "", opt[["steps"]])
  n <- unlist(lapply(strsplit(s, ",")[[1]], function(p) if (grepl("-", p)) { r <- as.integer(strsplit(p, "-")[[1]]); seq(r[1], r[2]) } else as.integer(p)))
  n <- sort(unique(n))
  if (any(is.na(n)) || any(n < 1 | n > length(COMPARE_STEPS))) stop("--steps must use numbers 1 to ", length(COMPARE_STEPS))
  n
})

# ---- Parameters -------------------------------------------------------------------
defaults <- load_compare_defaults()
params <- defaults
if (!is.null(opt[["based_on"]])) {
  base <- file.path(opt[["based_on"]], "params.yaml")
  if (!file.exists(base)) stop("--based_on folder has no params.yaml: ", opt[["based_on"]])
  params <- merge_compare_params(params, yaml::read_yaml(base))
  params$based_on <- basename(normalizePath(opt[["based_on"]]))
}
if (!is.null(opt[["config"]])) {
  cfg <- resolve_path(opt[["config"]])
  if (!file.exists(cfg)) stop("Config file not found: ", opt[["config"]])
  params <- merge_compare_params(params, yaml::read_yaml(cfg))
}
if (!is.null(opt[["studies"]])) {
  params$studies <- lapply(strsplit(opt[["studies"]], ",")[[1]], function(p) {
    kv <- strsplit(trimws(p), "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("--studies entries must be label=path, got '", p, "'")
    path <- resolve_path(kv[2])
    if (dir.exists(path)) list(label = kv[1], run_dir = path) else list(label = kv[1], de_file = path)
  })
}
if (!is.null(opt[["run_name"]])) params$run_name <- opt[["run_name"]]
params <- merge_compare_params(params, parse_set_string(opt[["set"]]))
for (k in COMPARE_LIST_PARAMS) if (k %in% names(params)) params[[k]] <- as_list_param(params[[k]])
unknown <- setdiff(names(params), c(names(defaults), "based_on"))
if (length(unknown)) log_warn("Unknown parameter(s) ignored (typo?): ", paste(unknown, collapse = ", "))

# ---- Run folder with a copy of each study's results -----------------------------------
runs_dir <- opt[["runs_dir"]] %||% file.path(REPO_ROOT, "runs")
dir.create(runs_dir, showWarnings = FALSE, recursive = TRUE)
run_dir <- file.path(normalizePath(runs_dir), paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "_", make.names(params$run_name %||% "compare")))
while (dir.exists(run_dir)) { Sys.sleep(1); run_dir <- paste0(run_dir, "b") }
dir.create(file.path(run_dir, "inputs"), recursive = TRUE)
dir.create(file.path(run_dir, "logs"))

# Copy what C01-C05 read, so the comparison can be rerun from this folder alone.
# A --based_on run already points at its own inputs/ copies.
copy_study <- function(s, i) {
  label <- as.character(s$label %||% paste0("study", i))
  dest <- file.path(run_dir, "inputs", sprintf("%02d_%s", i, safe_label(label)))
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  src_run <- resolve_path(s$run_dir %||% s$input_dir %||% "")
  if (!is.null(s$run_dir) || (!is.null(s$input_dir) && is.null(s$de_file))) {
    if (!dir.exists(src_run)) return(s)   # C01 reports the missing folder
    files <- c("04_differential_expression/de_results_all.csv", "04_differential_expression/de_summary.json",
               "01_validate_inputs/resolved_params.yaml", "provenance.json", "run_summary.json",
               list.files(src_run, pattern = "^gsea\\.csv$", recursive = TRUE))
    files <- unique(files[grepl("^(04_|01_|06_enrichment/|provenance|run_summary)", files)])
    for (f in files) {
      if (!file.exists(file.path(src_run, f))) next
      dir.create(dirname(file.path(dest, f)), recursive = TRUE, showWarnings = FALSE)
      file.copy(file.path(src_run, f), file.path(dest, f), overwrite = TRUE)
    }
    s$original_run_dir <- s$original_run_dir %||% normalizePath(src_run)
    s$run_dir <- NULL
  } else if (!is.null(s$de_file)) {
    src <- resolve_path(s$de_file)
    if (!file.exists(src)) return(s)
    target <- file.path(dest, paste0("de_file.", sub("^.*?\\.", "", basename(src))))
    file.copy(src, target, overwrite = TRUE)
    s$original_de_file <- s$original_de_file %||% normalizePath(src)
    s$de_file <- normalizePath(target)
  }
  s$input_dir <- normalizePath(dest)
  s
}
params$studies <- lapply(seq_along(params$studies), function(i) copy_study(params$studies[[i]], i))
write_yaml_file(params, file.path(run_dir, "params.yaml"),
                header = c("Parameters requested for this comparison (defaults + config + overrides).",
                           "Step C01 records how genes were matched and ranked: see C01_harmonise_studies/resolved_params.yaml"))
log_info("New comparison run: ", run_dir)

# ---- Steps --------------------------------------------------------------------------
rscript <- file.path(R.home("bin"), "Rscript")
records <- list()
overall <- "success"
for (i in steps) {
  st <- COMPARE_STEPS[[i]]
  log_file <- file.path(run_dir, "logs", paste0(st$name, ".log"))
  log_info(sprintf("Step %d/%d: %s ...", i, length(COMPARE_STEPS), st$name))
  t0 <- Sys.time()
  code <- system(paste(shQuote(rscript), shQuote(compare_step_script(st$name)), "--run_dir", shQuote(run_dir),
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
code_bank <- list(version = CODE_BANK_VERSION, module = "compare",
                  commit = if (length(commit) && is.null(attr(commit, "status"))) commit[1] else NULL,
                  uncommitted_changes = length(git("status", "--porcelain", "--untracked-files=no")) > 0)
resolved_path <- file.path(run_dir, "C01_harmonise_studies", "resolved_params.yaml")
all_pkgs <- list()
for (s in names(records)) {
  p <- read_json_file(file.path(run_dir, "provenance", paste0(s, ".json")))
  if (!is.null(p$packages)) all_pkgs <- modifyList(all_pkgs, as.list(p$packages))
}
write_json_file(list(
  run_id = basename(run_dir), module = "compare", status = overall, finished = timestamp_now(),
  code_bank = code_bank, r_version = R.version.string, platform = R.version$platform, packages = all_pkgs,
  parameters_requested = yaml::read_yaml(file.path(run_dir, "params.yaml")),
  parameters_resolved = if (file.exists(resolved_path)) yaml::read_yaml(resolved_path) else NULL,
  input_md5 = md5_of(list.files(file.path(run_dir, "inputs"), full.names = TRUE, recursive = TRUE)),
  steps = records,
  rerun_command = sprintf("Rscript compare/scripts/run_compare.R --based_on runs/%s", basename(run_dir))
), file.path(run_dir, "provenance.json"))

val <- read_json_file(file.path(run_dir, "C01_harmonise_studies", "validation_report.json"))
gs <- read_json_file(file.path(run_dir, "C02_compare_genes", "gene_comparison.json"))
rr <- read_json_file(file.path(run_dir, "C03_rrho", "rrho_summary.json"))
pw <- read_json_file(file.path(run_dir, "C04_compare_pathways", "pathway_comparison.json"))
warn <- if (!is.null(val) && length(val$checks)) {
  ch <- val$checks
  ch <- ch[ch$status %in% c("warning", "fail"), , drop = FALSE]
  if (NROW(ch)) as.list(paste0("[validation:", ch$id, "] ", ch$message)) else list()
} else list()
if (!is.null(gs$warnings)) warn <- c(warn, as.list(gs$warnings))
write_json_file(list(
  run_id = basename(run_dir), module = "compare", status = overall,
  steps = lapply(records, `[[`, "status"),
  studies = val$studies, match_method = val$match_method, rank_by = val$rank_by,
  significance_rule = gs$significance_rule,
  pairs = gs$pairs, rrho = rr$pairs, pathways = pw$pairs,
  warnings = warn
), file.path(run_dir, "run_summary.json"))
log_info("Comparison finished: ", overall, " -> ", run_dir)
quit(save = "no", status = if (overall == "failed") 1 else 0)
