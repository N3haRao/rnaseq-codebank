#!/usr/bin/env Rscript
# =============================================================================
# SU02_suppa_events_psi.R  --  events from the annotation, PSI per event
# =============================================================================
# 1. suppa.py generateEvents -f ioe -e <event_types>  (strict boundaries)
#    -> one .ioe file per event type, merged into events_all.ioe
# 2. suppa.py psiPerEvent for each group, from SU01's TPM tables
#
# Outputs (in <run_dir>/SU02_suppa_events_psi/): events_*_strict.ioe, events_all.ioe,
# <group>.psi, events_summary.json, suppa_*.log
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "suppa/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "suppa", "R", "suppa_utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "SU02_suppa_events_psi"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character"))))
params <- load_suppa_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
in_dir <- file.path(run_dir, "SU01_validate_suppa")
cmd <- strsplit(params$suppa_command, " ", fixed = TRUE)[[1]]

# ---- 1. Events ----------------------------------------------------------------------------
log_info("Finding alternative-splicing events in the annotation (", paste(params$event_types, collapse = ", "), ")")
prefix <- file.path(out_dir, "events")
code <- run_suppa(cmd, c("generateEvents", "-i", params$gtf_used, "-o", prefix, "-f", "ioe", "-e", params$event_types),
                  file.path(out_dir, "suppa_generateEvents.log"))
if (code != 0) fail("SUPPA generateEvents failed (see suppa_generateEvents.log).")
ioe <- sort(list.files(out_dir, pattern = "^events_.*_strict\\.ioe$", full.names = TRUE))
if (!length(ioe)) fail("generateEvents wrote no .ioe files.")
all_ioe <- file.path(out_dir, "events_all.ioe")
con <- file(all_ioe, "w")
for (i in seq_along(ioe)) {
  x <- readLines(ioe[i], warn = FALSE)
  writeLines(if (i == 1) x else x[-1], con)          # one header
}
close(con)
ev <- data.table::fread(all_ioe, data.table = FALSE)
if (!nrow(ev)) fail("No alternative-splicing events in this annotation (genes need two or more transcripts).")
types <- sub(":.*$", "", sub("^[^;]*;", "", ev$event_id))
by_type <- as.list(table(factor(types, levels = intersect(SUPPA_EVENT_ORDER, unique(types)))))
log_info(sprintf("%d events: %s", nrow(ev), paste(sprintf("%s %d", names(by_type), unlist(by_type)), collapse = ", ")))

# ---- 2. PSI per group -----------------------------------------------------------------------
groups <- c(reference = params$reference_level, test = params$test_level)
psi_files <- character(0)
for (k in names(groups)) {
  tpm_file <- file.path(in_dir, params$tpm_files[[k]])
  out <- file.path(out_dir, make.names(groups[[k]]))
  code <- run_suppa(cmd, c("psiPerEvent", "-i", all_ioe, "-e", tpm_file, "-o", out),
                    file.path(out_dir, sprintf("suppa_psiPerEvent_%s.log", make.names(groups[[k]]))))
  if (code != 0 || !file.exists(paste0(out, ".psi"))) fail("SUPPA psiPerEvent failed for ", groups[[k]], ".")
  psi_files[k] <- paste0(out, ".psi")
}
psi <- lapply(psi_files, read_suppa_matrix)
n_quant <- vapply(psi, function(m) sum(rowSums(is.na(m)) == 0), 0)
write_json_file(list(n_events = nrow(ev), events_by_type = by_type, event_types_requested = params$event_types,
                     events_with_psi_in_all_samples = as.list(setNames(n_quant, groups)),
                     psi_files = as.list(setNames(basename(psi_files), groups))),
                file.path(out_dir, "events_summary.json"))
write_step_provenance(run_dir, STEP, started, params_used = params[c("event_types")],
                      outputs = list.files(out_dir), notes = list(suppa_command = params$suppa_command))
log_info("Step SU02 done.")
