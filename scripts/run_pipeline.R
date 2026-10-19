#!/usr/bin/env Rscript
# =============================================================================
# run_pipeline.R  --  run the full analysis (or chosen steps) as one run
# =============================================================================
# This is the single entry point the website / agent calls. It:
#   1. creates a new run folder  runs/<YYYYmmdd-HHMMSS>_<run_name>/
#   2. merges parameters (defaults < --config < --based_on run < flags) and
#      saves them as params.yaml in the run folder
#   3. copies the input files into the run folder (so a run is self-contained)
#   4. runs each step script in its own R process, saving its log in logs/
#   5. stops at the first failed *critical* step (01-05); steps 06-07 are
#      allowed to fail without invalidating the core results
#   6. writes provenance.json (code commit, R, packages, parameters, input
#      checksums, step timings) and run_summary.json (headline results that
#      the chat agent reads first)
#
# Usage
#   # new run from a config file
#   Rscript scripts/run_pipeline.R --config config/examples/airway.yaml
#
#   # new run from flags only
#   Rscript scripts/run_pipeline.R --counts_file counts.csv --metadata_file samples.csv \
#       --condition genotype --reference_level WT --test_level KO --run_name ko_vs_wt
#
#   # rerun with a change: copy an earlier run's parameters, change one thing
#   Rscript scripts/run_pipeline.R --based_on runs/20261019-101500_airway --set "log2fc_cutoff=0.58"
#
#   # run only some steps of a new run (e.g. just validation)
#   Rscript scripts/run_pipeline.R --config my.yaml --steps 1
#
#   # resume an existing run from step 6 (same parameters)
#   Rscript scripts/run_pipeline.R --run_dir runs/20261019-101500_airway --steps 6-7
#
# Exit code: 0 = all critical steps passed, 1 = a critical step failed.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "scripts/x"), ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--config", type = "character", help = "Run config YAML (only the values that differ from config/default_params.yaml)"),
  make_option("--based_on", type = "character", help = "Start from the parameters of an earlier run folder"),
  make_option("--run_dir", type = "character", help = "Resume an EXISTING run folder (its saved parameters are used)"),
  make_option("--runs_dir", type = "character", default = NULL, help = "Where new run folders are created (default: <repo>/runs)"),
  make_option("--run_name", type = "character", help = "Short label for the run folder name"),
  make_option("--steps", type = "character", default = "1-7", help = "Steps to run, e.g. 1-7, 1, 4-5, 1,2,3 [default %default]"),
  make_option("--set", type = "character", help = "Overrides, e.g. \"padj_cutoff=0.01;covariates=[batch]\""),
  make_option("--counts_file", type = "character", help = "Counts table"),
  make_option("--metadata_file", type = "character", help = "Sample metadata table"),
  make_option("--condition", type = "character", help = "Metadata column with the groups"),
  make_option("--reference_level", type = "character", help = "Control / baseline group"),
  make_option("--test_level", type = "character", help = "Group compared to the reference"),
  make_option("--covariates", type = "character", help = "Comma-separated covariate columns"),
  make_option("--organism", type = "character", help = "auto | human | mouse | other")
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Run the RNA-seq code bank pipeline."))

# ---- Which steps --------------------------------------------------------------
parse_steps <- function(s) {
  nums <- unlist(lapply(strsplit(gsub(" ", "", s), ",")[[1]], function(p) {
    if (grepl("-", p)) { r <- as.integer(strsplit(p, "-")[[1]]); seq(r[1], r[2]) } else as.integer(p)
  }))
  nums <- sort(unique(nums))
  if (any(is.na(nums)) || any(nums < 1 | nums > length(PIPELINE_STEPS))) {
    stop("--steps must use numbers 1 to ", length(PIPELINE_STEPS), ", e.g. 1-7 or 4,5")
  }
  nums
}
step_nums <- parse_steps(opt$steps)

# ---- Parameters and run folder --------------------------------------------------
cli_keys <- c("counts_file", "metadata_file", "condition", "reference_level", "test_level", "covariates", "organism", "run_name")
cli <- opt[intersect(cli_keys, names(opt))]
cli <- cli[!vapply(cli, is.null, logical(1))]
overrides <- modifyList(cli, parse_set_string(opt$set))

if (!is.null(opt$run_dir)) {
  # Resume: keep the run's saved parameters exactly.
  run_dir <- normalizePath(opt$run_dir, mustWork = FALSE)
  if (!file.exists(file.path(run_dir, "params.yaml"))) stop("Not a run folder (no params.yaml): ", run_dir)
  if (length(overrides) || !is.null(opt$config)) {
    stop("--run_dir resumes a run with its saved parameters. To change parameters, start a new run with --based_on ", run_dir)
  }
  params <- yaml::read_yaml(file.path(run_dir, "params.yaml"))
  log_info("Resuming run ", run_dir)
} else {
  defaults <- load_defaults()
  params <- defaults
  if (!is.null(opt$based_on)) {
    base <- file.path(opt$based_on, "params.yaml")
    if (!file.exists(base)) stop("--based_on folder has no params.yaml: ", opt$based_on)
    params <- modifyList(params, yaml::read_yaml(base))
    params$based_on <- basename(normalizePath(opt$based_on))
    # point at the earlier run's input copies
    for (k in c("counts_file", "metadata_file")) {
      if (!is_blank(params[[k]]) && !file.exists(params[[k]])) {
        alt <- file.path(opt$based_on, "inputs", basename(params[[k]]))
        if (file.exists(alt)) params[[k]] <- normalizePath(alt)
      }
    }
  }
  if (!is.null(opt$config)) {
    cfg_path <- resolve_path(opt$config)
    if (!file.exists(cfg_path)) stop("Config file not found: ", opt$config)
    cfg <- yaml::read_yaml(cfg_path) %||% list()
    params <- modifyList(params, cfg)
  }
  params <- modifyList(params, overrides)
  params <- normalise_params(params)

  unknown <- setdiff(names(params), c(names(defaults), "based_on"))
  if (length(unknown)) {
    log_warn("Unknown parameter(s) ignored by the steps (typo?): ", paste(unknown, collapse = ", "))
  }

  runs_dir <- opt$runs_dir %||% file.path(REPO_ROOT, "runs")
  dir.create(runs_dir, showWarnings = FALSE, recursive = TRUE)
  run_id <- paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "_", make.names(params$run_name %||% "analysis"))
  run_dir <- file.path(normalizePath(runs_dir), run_id)
  while (dir.exists(run_dir)) { Sys.sleep(1); run_dir <- paste0(run_dir, "b") }
  dir.create(file.path(run_dir, "inputs"), recursive = TRUE)
  dir.create(file.path(run_dir, "logs"))

  # Copy inputs into the run so it is self-contained
  for (k in c("counts_file", "metadata_file")) {
    src <- resolve_path(params[[k]])
    if (!is_blank(src) && file.exists(src)) {
      dest <- file.path(run_dir, "inputs", basename(src))
      file.copy(src, dest, overwrite = TRUE)
      params[[paste0("original_", k)]] <- src
      params[[k]] <- normalizePath(dest)
    }
  }
  write_yaml_file(params, file.path(run_dir, "params.yaml"),
                  header = c("Parameters requested for this run (defaults + config + overrides).",
                             "Step 01 fills in anything left open: see 01_validate_inputs/resolved_params.yaml"))
  log_info("New run: ", run_dir)
}
dir.create(file.path(run_dir, "logs"), showWarnings = FALSE)

# ---- Code-bank version (git) ------------------------------------------------------
git_info <- function() {
  g <- function(...) tryCatch(suppressWarnings(system2("git", c("-C", shQuote(REPO_ROOT), ...),
                                                       stdout = TRUE, stderr = FALSE)),
                              error = function(e) character(0))
  commit <- g("rev-parse", "HEAD")
  if (!length(commit) || !is.null(attr(commit, "status"))) {
    return(list(commit = NULL, note = "Not a git repository (or git not installed): commit not recorded."))
  }
  status <- g("status", "--porcelain", "--untracked-files=no")
  list(commit = commit[1],
       branch = g("rev-parse", "--abbrev-ref", "HEAD")[1],
       uncommitted_changes = length(status) > 0,
       remote = g("config", "--get", "remote.origin.url")[1] %||% NULL,
       tag = { t <- g("describe", "--tags", "--exact-match"); if (length(t) && is.null(attr(t, "status"))) t[1] else NULL })
}
code_bank <- c(list(version = CODE_BANK_VERSION), git_info())
if (isTRUE(code_bank$uncommitted_changes)) {
  log_warn("The code bank has uncommitted changes; this run cannot be reproduced exactly from the commit alone.")
}

# ---- Run the steps ------------------------------------------------------------------
rscript <- file.path(R.home("bin"), "Rscript")
step_records <- list()
overall <- "success"
run_started <- Sys.time()
for (i in step_nums) {
  st <- PIPELINE_STEPS[[i]]
  script <- file.path(REPO_ROOT, "scripts", paste0(st$name, ".R"))
  log_file <- file.path(run_dir, "logs", paste0(st$name, ".log"))
  log_info(sprintf("Step %d/%d: %s ...", i, length(PIPELINE_STEPS), st$name))
  t0 <- Sys.time()
  cmd <- paste(shQuote(rscript), shQuote(script), "--run_dir", shQuote(run_dir), ">", shQuote(log_file), "2>&1")
  code <- system(cmd)
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  ok <- identical(as.integer(code), 0L)
  step_records[[st$name]] <- list(status = if (ok) "success" else "failed", exit_code = code,
                                  seconds = secs, critical = st$critical,
                                  log = file.path("logs", basename(log_file)))
  if (ok) {
    log_info(sprintf("  done in %.1f s", secs))
  } else {
    tail_lines <- tail(readLines(log_file, warn = FALSE), 15)
    log_msg(sprintf("  FAILED (exit code %s). Last lines of %s:", code, log_file), level = "ERROR")
    cat(paste0("    | ", tail_lines, "\n"), sep = "")
    if (st$critical) { overall <- "failed"; break }
    overall <- "success_with_errors"
  }
}

# ---- Provenance and run summary -------------------------------------------------------
resolved_path <- file.path(run_dir, "01_validate_inputs", "resolved_params.yaml")
resolved <- if (file.exists(resolved_path)) yaml::read_yaml(resolved_path) else NULL
step_prov <- lapply(setNames(nm = names(step_records)), function(s) read_json_file(file.path(run_dir, "provenance", paste0(s, ".json"))))
all_pkgs <- list()
for (p in step_prov) if (!is.null(p$packages)) all_pkgs <- modifyList(all_pkgs, as.list(p$packages))
prov <- list(
  run_id = basename(run_dir),
  status = overall,
  started = format(run_started, "%Y-%m-%d %H:%M:%S %Z"),
  finished = timestamp_now(),
  code_bank = code_bank,
  r_version = R.version.string,
  platform = R.version$platform,
  packages = all_pkgs,
  parameters_requested = yaml::read_yaml(file.path(run_dir, "params.yaml")),
  parameters_resolved = resolved,
  input_md5 = md5_of(list.files(file.path(run_dir, "inputs"), full.names = TRUE)),
  steps = step_records,
  rerun_command = sprintf("Rscript scripts/run_pipeline.R --based_on runs/%s", basename(run_dir))
)
write_json_file(prov, file.path(run_dir, "provenance.json"))

val <- read_json_file(file.path(run_dir, "01_validate_inputs", "validation_report.json"))
qc <- read_json_file(file.path(run_dir, "02_qc", "qc_summary.json"))
ex <- read_json_file(file.path(run_dir, "03_exploratory", "exploratory_summary.json"))
de <- read_json_file(file.path(run_dir, "04_differential_expression", "de_summary.json"))
en <- read_json_file(file.path(run_dir, "06_enrichment", "enrichment_summary.json"))
warnings_list <- c(
  if (!is.null(val) && length(val$checks)) {
    ch <- val$checks; ch <- ch[ch$status %in% c("warning", "fail"), , drop = FALSE]
    if (NROW(ch)) paste0("[validation:", ch$id, "] ", ch$message)
  },
  if (!is.null(qc) && NROW(qc$flags)) paste0("[qc:", qc$flags$flag, "] ", qc$flags$message),
  if (!is.null(ex) && NROW(ex$flagged_samples)) paste0("[pca:", ex$flagged_samples$flag, "] ", ex$flagged_samples$message)
)
summary <- list(
  run_id = basename(run_dir),
  status = overall,
  steps = lapply(step_records, `[[`, "status"),
  comparison = if (!is.null(resolved)) sprintf("%s vs %s (column '%s')", resolved$test_level, resolved$reference_level, resolved$condition),
  design_formula = resolved$design_formula,
  validation_status = val$status,
  warnings = as.list(warnings_list),
  n_samples = val$summary$n_samples,
  group_sizes = val$summary$group_sizes,
  n_genes_tested = de$n_genes_tested,
  significance_rule = de$significance_rule,
  n_significant = de$n_significant, n_up = de$n_up, n_down = de$n_down,
  top_up = de$top_up, top_down = de$top_down,
  pca_percent_variance = ex$percent_variance,
  enrichment = if (is.null(en)) NULL else if (!identical(en$status, "done")) list(status = en$status, reason = en$reason)
               else lapply(en$collections, function(cc) list(gsea_significant = cc$gsea$n_significant,
                                                             top_higher_in_test = cc$gsea$top_higher_in_test,
                                                             top_higher_in_reference = cc$gsea$top_higher_in_reference)),
  key_files = list(
    validation_report = "01_validate_inputs/validation_report.md",
    qc_metrics = "02_qc/qc_metrics.csv",
    pca_plot = "03_exploratory/pca_plot.png",
    de_results = "04_differential_expression/de_results_all.csv",
    de_significant = "04_differential_expression/de_results_significant.csv",
    volcano = "05_de_plots/volcano_plot.png",
    enrichment = "06_enrichment/enrichment_summary.json",
    report = "07_report/report.html",
    methods = "07_report/methods.md",
    provenance = "provenance.json"
  )
)
write_json_file(summary, file.path(run_dir, "run_summary.json"))

# ---- Console summary ----------------------------------------------------------------
cat("\n================================================================\n")
cat(sprintf(" Run:     %s\n Status:  %s\n", basename(run_dir), toupper(overall)))
if (!is.null(resolved)) cat(sprintf(" Compare: %s\n Design:  %s\n", summary$comparison, resolved$design_formula))
if (!is.null(de)) cat(sprintf(" DEGs:    %d significant (%d up, %d down); %s\n", de$n_significant, de$n_up, de$n_down, de$significance_rule))
if (length(warnings_list)) cat(sprintf(" Warnings: %d (see run_summary.json)\n", length(warnings_list)))
cat(sprintf(" Folder:  %s\n", run_dir))
if (file.exists(file.path(run_dir, "07_report", "report.html"))) cat(" Report:  07_report/report.html\n")
cat("================================================================\n")
quit(save = "no", status = if (overall == "failed") 1 else 0)
