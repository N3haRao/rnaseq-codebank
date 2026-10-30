#!/usr/bin/env Rscript
# =============================================================================
# S01_validate_rmats.R  --  check the rMATS output tables and standardise them
# =============================================================================
# Reads <EVENT>.MATS.<count_type>.txt for each event type, checks that they are
# complete rMATS statistical output, and writes one standard table of all events.
#
# Outputs (in <run_dir>/S01_validate_rmats/):
#   validation_report.json / .md   every check (pass / warning / fail) in plain language
#   events_validated.csv        one row per event, all event types
#   resolved_params.yaml           group names, replicate numbers, organism, event types found
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "splicing/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "splicing", "R", "splicing_utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "S01_validate_rmats"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "Splicing step S01: check rMATS tables."))
if (is_blank(opt$run_dir)) fail("--run_dir is required.")
params <- load_splicing_params(opt$run_dir, use_resolved = FALSE)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
log_info("Step S01: checking the rMATS tables")

checks <- list()
check <- function(id, status, message) {
  checks[[length(checks) + 1]] <<- list(id = id, status = status, message = message)
  (if (status == "pass") log_info else log_warn)(sprintf("[%s] %s: %s", status, id, message))
}
finish <- function() {
  statuses <- vapply(checks, `[[`, "", "status")
  overall <- if (any(statuses == "fail")) "fail" else if (any(statuses == "warning")) "warning" else "pass"
  rep <- list(status = overall, n_fail = sum(statuses == "fail"), n_warning = sum(statuses == "warning"), checks = checks)
  write_json_file(rep, file.path(out_dir, "validation_report.json"))
  md <- c("# rMATS input check", "", sprintf("Overall: **%s**", overall), "",
          vapply(checks, function(c) sprintf("- **%s** (%s): %s", c$id, c$status, c$message), ""))
  writeLines(md, file.path(out_dir, "validation_report.md"))
  overall
}

# ---- Files ---------------------------------------------------------------------
rmats_dir <- resolve_path(params$rmats_dir)
if (is.null(rmats_dir) || !dir.exists(rmats_dir)) {
  check("files_found", "fail", sprintf("The rMATS folder '%s' was not found. Set rmats_dir to the folder with SE.MATS.JCEC.txt and the other tables.", params$rmats_dir %||% ""))
  finish(); fail("No rMATS folder.")
}
count_type <- toupper(params$count_type %||% "JCEC")
if (!count_type %in% c("JC", "JCEC")) {
  check("count_type", "fail", "count_type must be JCEC or JC."); finish(); fail("Bad count_type.")
}
wanted <- intersect(params$event_types, EVENT_TYPES)
if (!length(wanted)) {
  check("event_types", "fail", paste("event_types must include some of:", paste(EVENT_TYPES, collapse = ", ")))
  finish(); fail("No event types.")
}
paths <- setNames(file.path(rmats_dir, sprintf("%s.MATS.%s.txt", wanted, count_type)), wanted)
found <- file.exists(paths)
if (!any(found)) {
  other <- list.files(rmats_dir, pattern = "\\.MATS\\.(JC|JCEC)\\.txt$")
  check("files_found", "fail", sprintf("No %s tables (e.g. SE.MATS.%s.txt) in %s.%s", count_type, count_type, rmats_dir,
        if (length(other)) paste0(" Found instead: ", paste(head(other, 5), collapse = ", "), ". Check count_type.") else
          " Only the 'fromGTF' files? Then rMATS was run without statistics (--statoff), which this module cannot use."))
  finish(); fail("No rMATS tables.")
}
if (all(found)) {
  check("files_found", "pass", sprintf("All %d event types found (%s, %s reads).", length(wanted), paste(wanted, collapse = ", "), count_type))
} else {
  check("files_found", "warning", sprintf("Missing: %s. The other event types are analysed.", paste(basename(paths[!found]), collapse = ", ")))
}

# ---- Read and check each table --------------------------------------------------
events <- list()
n_rep <- list()
for (et in names(paths)[found]) {
  df <- read_rmats_table(paths[[et]])
  need <- c(RMATS_COMMON, EVENT_COORDS[[et]])
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    check(paste0("columns_", et), "fail", sprintf("%s is missing columns: %s. Is it an rMATS 4 statistical output table?",
                                                 basename(paths[[et]]), paste(miss, collapse = ", ")))
    next
  }
  if (!nrow(df)) {
    check(paste0("events_", et), "warning", sprintf("%s has no events.", basename(paths[[et]])))
    next
  }
  if (all(is.na(suppressWarnings(as.numeric(df$FDR))))) {
    check(paste0("statistics_", et), "fail", sprintf("%s has no p-values / FDR. rMATS was probably run with --statoff.", basename(paths[[et]])))
    next
  }
  n_rep[[et]] <- c(max(n_values(df$IJC_SAMPLE_1)), max(n_values(df$IJC_SAMPLE_2)))
  ev <- standardise_events(df, et)
  bad <- is.na(ev$event_start) | is.na(ev$event_end) | ev$event_start > ev$event_end
  if (any(bad)) {
    check(paste0("coordinates_", et), "warning", sprintf("%d %s event(s) with unreadable coordinates were left out.", sum(bad), et))
    ev <- ev[!bad, ]
  }
  events[[et]] <- ev
}
if (!length(events)) { finish(); fail("No usable rMATS tables.") }
check("columns", "pass", sprintf("Tables readable: %s.", paste(names(events), collapse = ", ")))
all_events <- do.call(rbind, events)

# ---- Replicates and groups --------------------------------------------------------
reps <- do.call(rbind, n_rep)
r1 <- max(reps[, 1]); r2 <- max(reps[, 2])
if (r1 < 2 || r2 < 2) {
  check("replicates", "warning", sprintf("Group 1 has %d and group 2 has %d replicate(s). rMATS's test needs replicates to estimate variability; with fewer than 2 per group, treat results as exploratory.", r1, r2))
} else if (r1 < 3 || r2 < 3) {
  check("replicates", "warning", sprintf("Group 1 has %d and group 2 has %d replicates; at least 3 per group are recommended.", r1, r2))
} else {
  check("replicates", "pass", sprintf("Group 1: %d replicates, group 2: %d replicates.", r1, r2))
}
g1 <- params$group_1_name %||% "group_1"
g2 <- params$group_2_name %||% "group_2"
if (is_blank(params$group_1_name) || is_blank(params$group_2_name)) {
  check("group_names", "warning", "group_1_name / group_2_name are not set; results say group_1 (rMATS --b1) and group_2 (--b2). Set them to the real group names.")
} else {
  check("group_names", "pass", sprintf("Group 1 (rMATS --b1) = %s, group 2 (--b2) = %s. Positive dPSI = more inclusion in %s.", g1, g2, g1))
}
s1 <- params$group_1_samples; s2 <- params$group_2_samples
if (length(s1) && length(s1) != r1) { check("sample_names", "warning", sprintf("group_1_samples has %d names but there are %d replicates; generic names are used.", length(s1), r1)); s1 <- character(0) }
if (length(s2) && length(s2) != r2) { check("sample_names", "warning", sprintf("group_2_samples has %d names but there are %d replicates; generic names are used.", length(s2), r2)); s2 <- character(0) }
if (!length(s1)) s1 <- sprintf("%s_%d", g1, seq_len(r1))
if (!length(s2)) s2 <- sprintf("%s_%d", g2, seq_len(r2))

# ---- Organism (for gene sets) --------------------------------------------------------
organism <- params$organism %||% "auto"
if (organism == "auto") {
  ids <- head(all_events$gene_id[!is.na(all_events$gene_id)], 500)
  organism <- if (mean(grepl("^ENSMUSG", ids)) > 0.5) "mouse" else if (mean(grepl("^ENSG\\d", ids)) > 0.5) "human" else "other"
  check("organism", if (organism == "other") "warning" else "pass",
        sprintf("Organism from gene IDs: %s.%s", organism, if (organism == "other") " Pathway analysis with MSigDB is skipped; set organism or use gmt gene sets." else ""))
}
if (all(is.na(all_events$gene_symbol))) {
  check("gene_symbols", "warning", "No gene symbols in the rMATS tables (geneSymbol is empty); pathway analysis cannot match gene sets.")
}

status <- finish()
if (status == "fail") fail("The rMATS input check failed; see validation_report.md.")

data.table::fwrite(all_events, file.path(out_dir, "events_validated.csv"))
resolved <- list(rmats_dir = rmats_dir, count_type = count_type, event_types = names(events),
                 group_1_name = g1, group_2_name = g2, group_1_samples = s1, group_2_samples = s2,
                 n_replicates_1 = r1, n_replicates_2 = r2, organism = organism,
                 n_events = as.list(table(all_events$event_type)))
write_yaml_file(resolved, file.path(out_dir, "resolved_params.yaml"),
                header = "Values used by all later splicing steps (filled in by S01).")
write_step_provenance(run_dir, STEP, started, params_used = params[c("rmats_dir", "count_type", "event_types", "organism")],
                      inputs = unname(paths[found]), outputs = list.files(out_dir))
log_info(sprintf("Step S01 done: %d events (%s).", nrow(all_events),
                 paste(names(table(all_events$event_type)), table(all_events$event_type), sep = " ", collapse = ", ")))
