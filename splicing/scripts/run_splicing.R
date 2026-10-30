#!/usr/bin/env Rscript
# =============================================================================
# run_splicing.R  --  alternative-splicing analysis from rMATS output (one run)
# =============================================================================
# The entry point of the splicing module, built like scripts/run_pipeline.R:
#   1. creates a run folder  runs/<YYYYmmdd-HHMMSS>_<run_name>/
#   2. merges parameters (splicing defaults < --config < --based_on run < flags)
#      and saves them as params.yaml
#   3. copies the rMATS tables into inputs/ so the run is self-contained
#   4. runs S01-S04, each in its own R process with its own log
#   5. writes provenance.json and run_summary.json
#
# Usage
#   Rscript splicing/scripts/run_splicing.R --rmats_dir path/to/rmats_output \
#       --group_1_name KO --group_2_name WT --run_name ko_vs_wt_splicing
#   Rscript splicing/scripts/run_splicing.R --config my_splicing.yaml
#   Rscript splicing/scripts/run_splicing.R --based_on runs/<run> --set "delta_psi_cutoff=0.2"
#
# rMATS itself needs BAM files; later/scripts/U04_rmats.R runs it. This module
# starts from its output tables, so it works on results made anywhere.
# Exit code: 0 = critical steps passed, 1 = a critical step failed.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "splicing/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "splicing", "R", "splicing_utils.R"))
suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config", type = "character", help = "Run config YAML (values that differ from splicing/config/default_params.yaml)"),
  make_option("--based_on", type = "character", help = "Start from the parameters of an earlier splicing run folder"),
  make_option("--runs_dir", type = "character", help = "Where run folders are created (default: <repo>/runs)"),
  make_option("--run_name", type = "character", help = "Short label for the run folder name"),
  make_option("--steps", type = "character", default = "1-4", help = "Steps to run, e.g. 1-4 or 2-4 [default %default]"),
  make_option("--set", type = "character", help = "Overrides, e.g. \"delta_psi_cutoff=0.2;event_types=[SE, RI]\""),
  make_option("--rmats_dir", type = "character", help = "Folder with the rMATS tables (SE.MATS.JCEC.txt, ...)"),
  make_option("--group_1_name", type = "character", help = "rMATS sample 1 (--b1), the test group"),
  make_option("--group_2_name", type = "character", help = "rMATS sample 2 (--b2), the reference group")
), description = "Alternative-splicing analysis of rMATS output."))

steps <- local({
  s <- gsub(" ", "", opt$steps)
  n <- unlist(lapply(strsplit(s, ",")[[1]], function(p) if (grepl("-", p)) { r <- as.integer(strsplit(p, "-")[[1]]); seq(r[1], r[2]) } else as.integer(p)))
  n <- sort(unique(n))
  if (any(is.na(n)) || any(n < 1 | n > length(SPLICING_STEPS))) stop("--steps must use numbers 1 to ", length(SPLICING_STEPS))
  n
})

# ---- Parameters -------------------------------------------------------------------
defaults <- load_splicing_defaults()
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
cli <- opt[intersect(c("rmats_dir", "group_1_name", "group_2_name", "run_name"), names(opt))]
params <- modifyList(params, c(cli[!vapply(cli, is.null, logical(1))], parse_set_string(opt$set)))
for (k in SPLICING_LIST_PARAMS) if (k %in% names(params)) params[[k]] <- as_list_param(params[[k]])
unknown <- setdiff(names(params), c(names(defaults), "based_on", "original_rmats_dir"))
if (length(unknown)) log_warn("Unknown parameter(s) ignored (typo?): ", paste(unknown, collapse = ", "))

# ---- Run folder with a copy of the rMATS tables ----------------------------------------
runs_dir <- opt$runs_dir %||% file.path(REPO_ROOT, "runs")
dir.create(runs_dir, showWarnings = FALSE, recursive = TRUE)
run_dir <- file.path(normalizePath(runs_dir), paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "_", make.names(params$run_name %||% "splicing")))
while (dir.exists(run_dir)) { Sys.sleep(1); run_dir <- paste0(run_dir, "b") }
dir.create(file.path(run_dir, "inputs", "rmats"), recursive = TRUE)
dir.create(file.path(run_dir, "logs"))
src <- resolve_path(params$rmats_dir)
if (!is.null(src) && dir.exists(src)) {
  tables <- list.files(src, pattern = "^(SE|MXE|A3SS|A5SS|RI)\\.MATS\\.(JC|JCEC)\\.txt$", full.names = TRUE)
  file.copy(tables, file.path(run_dir, "inputs", "rmats"))
  params$original_rmats_dir <- params$original_rmats_dir %||% src
  params$rmats_dir <- normalizePath(file.path(run_dir, "inputs", "rmats"))
}
write_yaml_file(params, file.path(run_dir, "params.yaml"),
                header = c("Parameters requested for this splicing run (defaults + config + overrides).",
                           "Step S01 fills in group names, replicates and organism: see S01_validate_rmats/resolved_params.yaml"))
log_info("New splicing run: ", run_dir)

# ---- Steps --------------------------------------------------------------------------
rscript <- file.path(R.home("bin"), "Rscript")
records <- list()
overall <- "success"
for (i in steps) {
  st <- SPLICING_STEPS[[i]]
  log_file <- file.path(run_dir, "logs", paste0(st$name, ".log"))
  log_info(sprintf("Step %d/%d: %s ...", i, length(SPLICING_STEPS), st$name))
  t0 <- Sys.time()
  code <- system(paste(shQuote(rscript), shQuote(splicing_step_script(st$name)), "--run_dir", shQuote(run_dir),
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
code_bank <- list(version = CODE_BANK_VERSION, module = "splicing",
                  commit = if (length(commit) && is.null(attr(commit, "status"))) commit[1] else NULL,
                  uncommitted_changes = length(git("status", "--porcelain", "--untracked-files=no")) > 0)
resolved_path <- file.path(run_dir, "S01_validate_rmats", "resolved_params.yaml")
all_pkgs <- list()
for (s in names(records)) {
  p <- read_json_file(file.path(run_dir, "provenance", paste0(s, ".json")))
  if (!is.null(p$packages)) all_pkgs <- modifyList(all_pkgs, as.list(p$packages))
}
write_json_file(list(
  run_id = basename(run_dir), module = "splicing", status = overall, finished = timestamp_now(),
  code_bank = code_bank, r_version = R.version.string, platform = R.version$platform, packages = all_pkgs,
  parameters_requested = yaml::read_yaml(file.path(run_dir, "params.yaml")),
  parameters_resolved = if (file.exists(resolved_path)) yaml::read_yaml(resolved_path) else NULL,
  input_md5 = md5_of(list.files(file.path(run_dir, "inputs", "rmats"), full.names = TRUE)),
  steps = records,
  rerun_command = sprintf("Rscript splicing/scripts/run_splicing.R --based_on runs/%s", basename(run_dir))
), file.path(run_dir, "provenance.json"))

val <- read_json_file(file.path(run_dir, "S01_validate_rmats", "validation_report.json"))
ss <- read_json_file(file.path(run_dir, "S02_splicing_results", "splicing_summary.json"))
warn <- if (!is.null(val) && length(val$checks)) {
  ch <- val$checks
  ch <- ch[ch$status %in% c("warning", "fail"), , drop = FALSE]
  if (NROW(ch)) as.list(paste0("[validation:", ch$id, "] ", ch$message)) else list()
} else list()
write_json_file(list(
  run_id = basename(run_dir), module = "splicing", status = overall,
  steps = lapply(records, `[[`, "status"),
  comparison = ss$comparison, significance_rule = ss$significance_rule, sign_convention = ss$sign_convention,
  n_events_tested = ss$n_events_tested, n_significant = ss$n_significant, n_up = ss$n_up, n_down = ss$n_down,
  n_genes_significant = ss$n_genes_significant, by_event_type = ss$by_event_type, top_events = ss$top_events,
  warnings = warn
), file.path(run_dir, "run_summary.json"))
log_info("Splicing run finished: ", overall, " -> ", run_dir)
quit(save = "no", status = if (overall == "failed") 1 else 0)
