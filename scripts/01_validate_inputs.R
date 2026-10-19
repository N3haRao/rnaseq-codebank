#!/usr/bin/env Rscript
# =============================================================================
# 01_validate_inputs.R  --  check the uploaded files before any analysis
# =============================================================================
# What it does
#   Reads the counts table and the sample metadata, checks them, explains any
#   problem in plain language, and writes clean, aligned copies for later steps.
#   It also resolves anything left open in the parameters (reference group,
#   organism, gene ID type, filtering threshold) and records those decisions.
#
# Checks (each gets an id, a status of pass / warning / fail, and a message)
#   files_found, counts_readable, gene_id_column, numeric_columns, excluded_samples,
#   missing_values, negative_values, count_values (raw counts vs TPM/log),
#   ensembl_versions, duplicate_gene_ids, all_zero_genes, library_sizes,
#   metadata_readable, sample_id_column, duplicate_samples, sample_matching,
#   condition_column, condition_values, reference_level, test_level,
#   replicates, covariates, design_full_rank, residual_df, level_names,
#   gene_ids_organism
#
# Outputs (in <run_dir>/01_validate_inputs/)
#   validation_report.json   every check, machine-readable (the agent reads this)
#   validation_report.md     the same, for people
#   counts_validated.csv     integer counts, gene IDs cleaned, samples in metadata order
#   metadata_validated.csv   one row per sample, first column renamed sample_id
#   column_types.json        which metadata columns are numeric
#   resolved_params.yaml     final parameters used by every later step
#
# Exit code: 0 when there are no failed checks, 1 otherwise.
#
# Usage
#   Rscript scripts/01_validate_inputs.R --run_dir runs/my_run
#   Rscript scripts/01_validate_inputs.R --run_dir runs/my_run \
#       --counts_file my_counts.csv --metadata_file my_samples.csv \
#       --condition genotype --reference_level WT --test_level KO
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "scripts/x"), ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "01_validate_inputs"
started <- Sys.time()

option_list <- list(
  make_option("--run_dir", type = "character", help = "Run folder (required)"),
  make_option("--config", type = "character", help = "Parameter file (default: <run_dir>/params.yaml)"),
  make_option("--set", type = "character", help = "Extra overrides, e.g. \"padj_cutoff=0.01;covariates=[batch]\""),
  make_option("--counts_file", type = "character", help = "Counts table (genes x samples)"),
  make_option("--metadata_file", type = "character", help = "Sample metadata table"),
  make_option("--gene_id_column", type = "character", help = "Gene ID column in counts (default: first)"),
  make_option("--sample_id_column", type = "character", help = "Sample ID column in metadata (default: first)"),
  make_option("--condition", type = "character", help = "Metadata column with the groups"),
  make_option("--reference_level", type = "character", help = "Control / baseline group"),
  make_option("--test_level", type = "character", help = "Group compared to the reference"),
  make_option("--covariates", type = "character", help = "Comma-separated covariate columns"),
  make_option("--exclude_samples", type = "character", help = "Comma-separated sample IDs to leave out"),
  make_option("--organism", type = "character", help = "auto | human | mouse | other"),
  make_option("--duplicate_gene_ids", type = "character", help = "sum | first | error"),
  make_option("--round_counts", type = "logical", help = "Round non-integer counts (TRUE/FALSE)"),
  make_option("--min_count", type = "integer", help = "Low-count filter: minimum reads"),
  make_option("--min_samples", type = "integer", help = "Low-count filter: minimum samples")
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Step 01: validate the counts and metadata files."))
params <- load_params(opt, use_resolved = FALSE)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
log_info("Step 01: validating inputs in ", run_dir)

# -----------------------------------------------------------------------------
# Check bookkeeping
# -----------------------------------------------------------------------------
checks <- list()
add_check <- function(id, status, message, details = NULL) {
  checks[[length(checks) + 1]] <<- list(id = id, status = status, message = message,
                                         details = details)
  lvl <- switch(status, pass = "INFO", warning = "WARN", fail = "ERROR")
  log_msg(sprintf("[%s] %s", toupper(status), message), level = lvl)
}
summary_info <- list()

show_some <- function(x, n = 5) {
  x <- unique(x)
  if (length(x) <= n) paste(x, collapse = ", ")
  else paste0(paste(head(x, n), collapse = ", "), ", ... (", length(x), " in total)")
}

finish <- function() {
  statuses <- vapply(checks, `[[`, character(1), "status")
  overall <- if (any(statuses == "fail")) "fail" else if (any(statuses == "warning")) "warning" else "pass"
  report <- list(step = STEP, status = overall, generated = timestamp_now(),
                 n_fail = sum(statuses == "fail"), n_warning = sum(statuses == "warning"),
                 summary = summary_info, checks = checks)
  write_json_file(report, file.path(out_dir, "validation_report.json"))

  icon <- c(pass = "PASS", warning = "WARNING", fail = "FAIL")
  md <- c("# Input validation report", "",
          sprintf("**Overall: %s**  (%d failed, %d warnings)", toupper(overall),
                  report$n_fail, report$n_warning), "",
          "| Status | Check | Message |", "| --- | --- | --- |",
          vapply(checks, function(ch) sprintf("| %s | %s | %s |", icon[[ch$status]], ch$id,
                                              gsub("\\|", "/", ch$message)), character(1)))
  writeLines(md, file.path(out_dir, "validation_report.md"))

  write_step_provenance(run_dir, STEP, started,
                        params_used = params[c("counts_file", "metadata_file", "gene_id_column",
                                               "sample_id_column", "condition", "reference_level",
                                               "test_level", "covariates", "organism",
                                               "duplicate_gene_ids", "round_counts", "exclude_samples")],
                        inputs = c(params$counts_file, params$metadata_file),
                        outputs = list.files(out_dir),
                        notes = list(overall_status = overall))
  if (overall == "fail") {
    log_msg("Validation FAILED. See ", file.path(out_dir, "validation_report.md"), level = "ERROR")
    quit(save = "no", status = 1)
  }
  log_info("Validation finished with status: ", toupper(overall))
  quit(save = "no", status = 0)
}

# -----------------------------------------------------------------------------
# 1. Files present
# -----------------------------------------------------------------------------
params$counts_file   <- resolve_path(params$counts_file)
params$metadata_file <- resolve_path(params$metadata_file)
missing_files <- c(
  if (is_blank(params$counts_file)) "counts_file is not set"
  else if (!file.exists(params$counts_file)) paste("counts file not found:", params$counts_file),
  if (is_blank(params$metadata_file)) "metadata_file is not set"
  else if (!file.exists(params$metadata_file)) paste("metadata file not found:", params$metadata_file)
)
if (length(missing_files)) {
  add_check("files_found", "fail", paste("Input files missing:", paste(missing_files, collapse = "; ")))
  finish()
}
add_check("files_found", "pass", "Both input files were found.")

# -----------------------------------------------------------------------------
# 2. Counts table
# -----------------------------------------------------------------------------
counts_df <- tryCatch(read_table_file(params$counts_file), error = function(e) e)
if (inherits(counts_df, "error") || ncol(counts_df) < 2 || nrow(counts_df) < 1) {
  add_check("counts_readable", "fail",
            paste0("The counts file could not be read as a table with gene rows and sample columns",
                   if (inherits(counts_df, "error")) paste0(" (", conditionMessage(counts_df), ")"), "."))
  finish()
}
add_check("counts_readable", "pass",
          sprintf("Counts file read: %d rows and %d columns.", nrow(counts_df), ncol(counts_df)))

# Gene ID column
gene_col <- params$gene_id_column %||% names(counts_df)[1]
if (!gene_col %in% names(counts_df)) {
  add_check("gene_id_column", "fail",
            sprintf("Gene ID column '%s' is not in the counts file. Columns are: %s.",
                    gene_col, show_some(names(counts_df), 10)))
  finish()
}
gene_ids <- trimws(as.character(counts_df[[gene_col]]))
if (is.numeric(counts_df[[gene_col]])) {
  add_check("gene_id_column", "fail",
            sprintf(paste("The gene ID column ('%s') contains numbers, not gene IDs. The first column",
                          "of the counts file should hold gene IDs (e.g. ENSG00000141510 or TP53)."),
                    gene_col))
  finish()
}
empty_id <- is.na(gene_ids) | gene_ids == ""
if (any(empty_id)) {
  add_check("gene_id_column", "warning",
            sprintf("%d rows have no gene ID and were removed.", sum(empty_id)))
  counts_df <- counts_df[!empty_id, , drop = FALSE]
  gene_ids <- gene_ids[!empty_id]
} else {
  add_check("gene_id_column", "pass", sprintf("Gene IDs taken from column '%s'.", gene_col))
}

# Sample (numeric) columns
other_cols <- setdiff(names(counts_df), gene_col)
is_num <- vapply(counts_df[other_cols], is.numeric, logical(1))
non_numeric <- other_cols[!is_num]
sample_cols <- other_cols[is_num]
if (length(sample_cols) == 0) {
  add_check("numeric_columns", "fail",
            "No columns of numbers were found in the counts file. Each sample should be one column of counts.")
  finish()
}
if (length(non_numeric)) {
  add_check("numeric_columns", "warning",
            sprintf(paste("These columns are not numbers and were ignored (they look like annotation",
                          "columns, not samples): %s."), show_some(non_numeric)))
} else {
  add_check("numeric_columns", "pass", sprintf("%d sample columns found.", length(sample_cols)))
}
counts <- as.matrix(counts_df[, sample_cols, drop = FALSE])
storage.mode(counts) <- "double"
colnames(counts) <- trimws(colnames(counts))

# Samples the user asked to leave out
excl <- as_list_param(params$exclude_samples)
if (length(excl)) {
  unknown <- setdiff(excl, colnames(counts))
  if (length(unknown)) {
    add_check("excluded_samples", "fail",
              sprintf("These samples were listed in exclude_samples but are not columns of the counts file: %s.",
                      paste(unknown, collapse = ", ")))
    finish()
  }
  counts <- counts[, !colnames(counts) %in% excl, drop = FALSE]
  add_check("excluded_samples", "warning",
            sprintf("Left out of the analysis at your request: %s.", paste(excl, collapse = ", ")))
}

# Missing values
n_na <- sum(is.na(counts))
if (n_na > 0) {
  bad_genes <- gene_ids[rowSums(is.na(counts)) > 0]
  add_check("missing_values", "fail",
            sprintf("The counts table has %d empty or non-numeric cells (genes: %s). Counts must be complete; use 0 for genes with no reads.",
                    n_na, show_some(bad_genes)))
  finish()
}
add_check("missing_values", "pass", "No missing values in the counts.")

# Negative values
if (any(counts < 0)) {
  add_check("negative_values", "fail",
            paste("The counts table contains negative numbers. Raw counts are never negative;",
                  "this usually means the data were log-transformed or batch-corrected.",
                  "Please upload the raw counts instead."))
  finish()
}
add_check("negative_values", "pass", "No negative values.")

# Raw counts vs normalised values
non_integer <- any(abs(counts - round(counts)) > 1e-6)
col_sums <- colSums(counts)
if (non_integer) {
  looks_tpm <- all(abs(col_sums - 1e6) / 1e6 < 0.02)
  looks_log <- max(counts) < 50
  if (looks_tpm) {
    add_check("count_values", "fail",
              paste("Every sample adds up to about 1,000,000, so these look like TPM or CPM values,",
                    "not raw counts. DESeq2 needs raw read counts (for example from featureCounts,",
                    "HTSeq, or the 'NumReads'/'expected_count' column of Salmon/RSEM)."))
    finish()
  }
  if (looks_log) {
    add_check("count_values", "fail",
              paste("All values are below 50 and have decimals, so they look log-transformed",
                    "(e.g. log2 CPM). Please upload the raw read counts instead."))
    finish()
  }
  if (isTRUE(params$round_counts)) {
    add_check("count_values", "warning",
              paste("Some counts have decimals. This is normal for estimated counts from Salmon or",
                    "RSEM; they were rounded to whole numbers, as DESeq2 requires."))
    counts <- round(counts)
  } else {
    add_check("count_values", "fail",
              "Some counts have decimals and round_counts is false. DESeq2 needs whole-number counts.")
    finish()
  }
} else {
  add_check("count_values", "pass", "Values are whole numbers, as expected for raw counts.")
}

# Ensembl version suffixes (ENSG00000141510.17 -> ENSG00000141510)
stripped <- strip_ensembl_version(gene_ids)
if (any(stripped != gene_ids)) {
  add_check("ensembl_versions", "pass",
            sprintf("Removed version numbers from %d Ensembl IDs (e.g. %s -> %s).",
                    sum(stripped != gene_ids), gene_ids[which(stripped != gene_ids)[1]],
                    stripped[which(stripped != gene_ids)[1]]))
  gene_ids <- stripped
}

# Duplicate gene IDs
dup <- duplicated(gene_ids)
if (any(dup)) {
  dup_ids <- unique(gene_ids[dup])
  how <- params$duplicate_gene_ids %||% "sum"
  if (how == "error") {
    add_check("duplicate_gene_ids", "fail",
              sprintf("%d gene IDs appear more than once (e.g. %s). Each gene should appear once.",
                      length(dup_ids), show_some(dup_ids)))
    finish()
  } else if (how == "first") {
    counts <- counts[!dup, , drop = FALSE]
    gene_ids <- gene_ids[!dup]
    add_check("duplicate_gene_ids", "warning",
              sprintf("%d gene IDs appear more than once (e.g. %s); the first row of each was kept.",
                      length(dup_ids), show_some(dup_ids)))
  } else {
    counts <- rowsum(counts, group = gene_ids, reorder = FALSE)
    gene_ids <- rownames(counts)
    add_check("duplicate_gene_ids", "warning",
              sprintf("%d gene IDs appear more than once (e.g. %s); their counts were added together.",
                      length(dup_ids), show_some(dup_ids)))
  }
} else {
  add_check("duplicate_gene_ids", "pass", "Every gene ID is unique.")
}
rownames(counts) <- gene_ids

n_zero <- sum(rowSums(counts) == 0)
add_check("all_zero_genes", "pass",
          sprintf("%d of %d genes have zero reads in every sample (they are removed by the low-count filter later).",
                  n_zero, nrow(counts)))

# Library sizes
lib <- colSums(counts)
min_lib <- params$min_library_size %||% 1e6
if (any(lib == 0)) {
  add_check("library_sizes", "fail",
            sprintf("These samples have no reads at all: %s. Remove them or check the file.",
                    show_some(names(lib)[lib == 0])))
  finish()
}
small <- names(lib)[lib < min_lib]
if (length(small)) {
  add_check("library_sizes", "warning",
            sprintf("These samples have fewer than %s reads, which is low for RNA-seq: %s.",
                    format(min_lib, big.mark = ",", scientific = FALSE),
                    paste(sprintf("%s (%s)", small, format(lib[small], big.mark = ",")), collapse = ", ")))
} else {
  add_check("library_sizes", "pass",
            sprintf("Library sizes range from %s to %s reads.",
                    format(min(lib), big.mark = ","), format(max(lib), big.mark = ",")))
}

# -----------------------------------------------------------------------------
# 3. Metadata table
# -----------------------------------------------------------------------------
meta <- tryCatch(read_table_file(params$metadata_file), error = function(e) e)
if (inherits(meta, "error") || nrow(meta) < 1 || ncol(meta) < 2) {
  add_check("metadata_readable", "fail",
            "The metadata file could not be read as a table with one row per sample and at least two columns (sample ID and group).")
  finish()
}
add_check("metadata_readable", "pass",
          sprintf("Metadata read: %d samples and %d columns (%s).", nrow(meta), ncol(meta),
                  show_some(names(meta), 8)))

sample_col <- params$sample_id_column %||% names(meta)[1]
if (!sample_col %in% names(meta)) {
  add_check("sample_id_column", "fail",
            sprintf("Sample ID column '%s' is not in the metadata. Columns are: %s.",
                    sample_col, show_some(names(meta), 10)))
  finish()
}
add_check("sample_id_column", "pass", sprintf("Sample IDs taken from metadata column '%s'.", sample_col))
meta[[sample_col]] <- trimws(as.character(meta[[sample_col]]))
meta <- meta[!(is.na(meta[[sample_col]]) | meta[[sample_col]] == ""), , drop = FALSE]
if (length(excl)) meta <- meta[!meta[[sample_col]] %in% excl, , drop = FALSE]

if (any(duplicated(meta[[sample_col]]))) {
  add_check("duplicate_samples", "fail",
            sprintf("These sample IDs appear more than once in the metadata: %s.",
                    show_some(meta[[sample_col]][duplicated(meta[[sample_col]])])))
  finish()
}
add_check("duplicate_samples", "pass", "Every sample ID in the metadata is unique.")

# Matching samples between the two files
in_counts_only <- setdiff(colnames(counts), meta[[sample_col]])
in_meta_only <- setdiff(meta[[sample_col]], colnames(counts))
if (length(in_counts_only)) {
  hint <- vapply(in_counts_only, function(s) {
    if (!length(in_meta_only)) return("")
    d <- adist(tolower(s), tolower(in_meta_only))[1, ]
    best <- in_meta_only[which.min(d)]
    if (min(d) <= max(3, nchar(s) %/% 3)) sprintf(" (closest metadata ID: '%s')", best) else ""
  }, character(1))
  add_check("sample_matching", "fail",
            paste0("These sample columns in the counts file have no row in the metadata: ",
                   paste0("'", in_counts_only, "'", hint, collapse = ", "),
                   ". Sample names must match exactly between the two files (including upper/lower case)."),
            details = list(in_counts_only = in_counts_only, in_metadata_only = in_meta_only))
  finish()
}
if (length(in_meta_only)) {
  add_check("sample_matching", "warning",
            sprintf("These metadata samples are not in the counts file and were ignored: %s.",
                    show_some(in_meta_only)))
} else {
  add_check("sample_matching", "pass",
            sprintf("All %d samples match between the counts and metadata files.", ncol(counts)))
}
meta <- meta[match(colnames(counts), meta[[sample_col]]), , drop = FALSE]
names(meta)[names(meta) == sample_col] <- "sample_id"
rownames(meta) <- meta$sample_id

# -----------------------------------------------------------------------------
# 4. Experimental design
# -----------------------------------------------------------------------------
cond <- params$condition
if (is_blank(cond) || !cond %in% names(meta)) {
  add_check("condition_column", "fail",
            sprintf("The group column '%s' is not in the metadata. Available columns: %s. Set 'condition' to the column that holds the groups you want to compare.",
                    cond %||% "(not set)", show_some(setdiff(names(meta), "sample_id"), 10)))
  finish()
}
add_check("condition_column", "pass", sprintf("Groups taken from metadata column '%s'.", cond))
cond_values <- trimws(as.character(meta[[cond]]))
if (any(is.na(cond_values) | cond_values == "")) {
  add_check("condition_values", "fail",
            sprintf("These samples have no value in '%s': %s.", cond,
                    show_some(meta$sample_id[is.na(cond_values) | cond_values == ""])))
  finish()
}
meta[[cond]] <- cond_values
group_sizes <- table(cond_values)
if (length(group_sizes) < 2) {
  add_check("condition_values", "fail",
            sprintf("Column '%s' has only one group (%s). At least two groups are needed to compare.",
                    cond, names(group_sizes)))
  finish()
}
add_check("condition_values", "pass",
          sprintf("Groups found: %s.", paste(sprintf("%s (n=%d)", names(group_sizes), group_sizes),
                                            collapse = ", ")))

levels_all <- names(group_sizes)
control_words <- c("control", "ctrl", "ctl", "wt", "wildtype", "wild_type", "wild-type", "untreated",
                   "untrt", "vehicle", "veh", "mock", "dmso", "pbs", "sham", "baseline", "normal",
                   "healthy", "naive", "scramble", "scr", "sicontrol", "shctrl", "nt", "neg")
ref <- params$reference_level
if (!is_blank(ref)) {
  if (!ref %in% levels_all) {
    add_check("reference_level", "fail",
              sprintf("Reference group '%s' is not in column '%s'. The groups are: %s.",
                      ref, cond, paste(levels_all, collapse = ", ")))
    finish()
  }
  add_check("reference_level", "pass", sprintf("Reference (control) group: '%s'.", ref))
} else {
  hit <- levels_all[tolower(levels_all) %in% control_words]
  ref <- if (length(hit) == 1) hit else sort(levels_all)[1]
  add_check("reference_level", "warning",
            sprintf("No reference group was given, so '%s' was chosen%s. Please confirm this is the control group.",
                    ref, if (length(hit) == 1) " because its name looks like a control" else " (first alphabetically)"))
}
test <- params$test_level
if (!is_blank(test)) {
  if (!test %in% levels_all || test == ref) {
    add_check("test_level", "fail",
              sprintf("Test group '%s' must be one of the groups in '%s' and differ from the reference. The groups are: %s.",
                      test, cond, paste(levels_all, collapse = ", ")))
    finish()
  }
  add_check("test_level", "pass", sprintf("Comparison: '%s' vs '%s'.", test, ref))
} else if (length(levels_all) == 2) {
  test <- setdiff(levels_all, ref)
  add_check("test_level", "pass", sprintf("Comparison: '%s' vs '%s'.", test, ref))
} else {
  add_check("test_level", "fail",
            sprintf("Column '%s' has %d groups (%s). Choose which group to compare with '%s' (test_level).",
                    cond, length(levels_all), paste(levels_all, collapse = ", "), ref))
  finish()
}

n_ref <- group_sizes[[ref]]; n_test <- group_sizes[[test]]
min_rep_warn <- params$min_replicates_warn %||% 3
if (min(n_ref, n_test) < 2) {
  add_check("replicates", "fail",
            sprintf("'%s' has %d sample(s) and '%s' has %d. Each compared group needs at least 2 biological replicates for a statistical test.",
                    ref, n_ref, test, n_test))
  finish()
} else if (min(n_ref, n_test) < min_rep_warn) {
  add_check("replicates", "warning",
            sprintf("'%s' has %d and '%s' has %d replicates. Fewer than %d per group gives low power: only large changes will be detected, and results are more sensitive to a single unusual sample.",
                    ref, n_ref, test, n_test, min_rep_warn))
} else {
  add_check("replicates", "pass", sprintf("Replicates: %d (%s) and %d (%s).", n_ref, ref, n_test, test))
}

# Covariates
covs <- as_list_param(params$covariates)
col_types <- list()
for (cv in covs) {
  if (!cv %in% names(meta)) {
    add_check("covariates", "fail",
              sprintf("Covariate '%s' is not in the metadata. Available columns: %s.", cv,
                      show_some(setdiff(names(meta), "sample_id"), 10)))
    finish()
  }
  if (cv == cond) {
    add_check("covariates", "fail", sprintf("'%s' is the group column and cannot also be a covariate.", cv))
    finish()
  }
  v <- meta[[cv]]
  if (any(is.na(v) | trimws(as.character(v)) == "")) {
    add_check("covariates", "fail",
              sprintf("Covariate '%s' is missing for: %s. Every sample needs a value.", cv,
                      show_some(meta$sample_id[is.na(v) | trimws(as.character(v)) == ""])))
    finish()
  }
  if (length(unique(v)) < 2) {
    add_check("covariates", "fail",
              sprintf("Covariate '%s' has the same value for every sample, so it cannot be adjusted for. Remove it from the covariates.", cv))
    finish()
  }
  if (is.numeric(v) && length(unique(v)) > 5) col_types[[cv]] <- "numeric"
}
if (length(covs)) {
  add_check("covariates", "pass",
            sprintf("Adjusting for: %s%s.", paste(covs, collapse = ", "),
                    if (length(col_types)) paste0(" (treated as numeric: ", paste(names(col_types), collapse = ", "), ")") else ""))
}

# Design: can every effect be estimated separately?
design_formula <- paste("~", paste(c(sprintf("`%s`", covs), sprintf("`%s`", cond)), collapse = " + "))
mm_data <- meta[, c(covs, cond), drop = FALSE]
for (cv in setdiff(c(covs, cond), names(col_types))) mm_data[[cv]] <- factor(as.character(mm_data[[cv]]))
mm <- model.matrix(as.formula(design_formula), data = mm_data)
rank <- qr(mm)$rank
if (rank < ncol(mm)) {
  conf <- vapply(covs, function(cv) {
    tab <- table(mm_data[[cv]], mm_data[[cond]])
    all(rowSums(tab > 0) == 1)
  }, logical(1))
  why <- if (any(conf)) {
    sprintf(" Each value of '%s' contains samples from only one group, so its effect cannot be separated from the group difference (they are confounded).",
            paste(covs[conf], collapse = "', '"))
  } else " Some combination of covariates and groups is redundant."
  add_check("design_full_rank", "fail",
            paste0("The design ", gsub("`", "", design_formula), " cannot be estimated.", why,
                   " Remove the confounded covariate, or note that this experiment cannot separate those effects."))
  finish()
}
add_check("design_full_rank", "pass", sprintf("Design %s can be estimated.", gsub("`", "", design_formula)))

resid_df <- nrow(meta) - rank
if (resid_df < 1) {
  add_check("residual_df", "fail",
            sprintf("There are %d samples but the design needs %d parameters, leaving nothing to estimate variability. Use fewer covariates.",
                    nrow(meta), rank))
  finish()
} else if (resid_df < 3) {
  add_check("residual_df", "warning",
            sprintf("Only %d degree(s) of freedom are left after fitting the design, so variability estimates rely heavily on sharing information across genes.", resid_df))
} else {
  add_check("residual_df", "pass", sprintf("%d residual degrees of freedom.", resid_df))
}

odd_names <- c(levels_all, unlist(lapply(setdiff(covs, names(col_types)), function(cv) unique(as.character(meta[[cv]])))))
odd_names <- odd_names[make.names(odd_names) != odd_names]
if (length(odd_names)) {
  add_check("level_names", "pass",
            sprintf("Some group names contain spaces or symbols (%s). This is fine: they are converted internally and reported with their original names.",
                    show_some(odd_names)))
}

# -----------------------------------------------------------------------------
# 5. Gene IDs and organism
# -----------------------------------------------------------------------------
det <- detect_gene_ids(gene_ids)
organism <- params$organism %||% "auto"
if (organism == "auto") {
  organism <- det$organism_guess
  add_check("gene_ids_organism", if (organism == "other") "warning" else "pass",
            sprintf("Gene IDs look like %s IDs from %s.%s", det$id_type,
                    if (organism == "other") "an organism other than human or mouse" else organism,
                    if (organism == "other") " Gene names and pathway analysis are available for human and mouse only; differential expression still runs." else ""))
} else {
  add_check("gene_ids_organism", "pass",
            sprintf("Organism set to '%s'; gene IDs look like %s IDs.", organism, det$id_type))
}

# -----------------------------------------------------------------------------
# 6. Write outputs
# -----------------------------------------------------------------------------
min_samples <- params$min_samples %||% min(n_ref, n_test)
resolved <- list(
  counts_file = params$counts_file,
  metadata_file = params$metadata_file,
  condition = cond,
  reference_level = ref,
  test_level = test,
  covariates = covs,
  exclude_samples = excl,
  numeric_covariates = names(col_types),
  organism = organism,
  gene_id_type = det$id_type,
  min_count = params$min_count %||% 10,
  min_samples = as.integer(min_samples),
  design_formula = gsub("`", "", design_formula)
)
write_yaml_file(resolved, file.path(out_dir, "resolved_params.yaml"),
                header = c("Written by step 01. Final values used by every later step.",
                           "Change them by editing the run config and re-running step 01."))

counts_out <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
write_table(counts_out, file.path(out_dir, "counts_validated.csv"))
write_table(meta, file.path(out_dir, "metadata_validated.csv"))
write_json_file(col_types, file.path(out_dir, "column_types.json"))

summary_info <- list(
  n_genes = nrow(counts), n_samples = ncol(counts),
  condition = cond, reference_level = ref, test_level = test,
  group_sizes = as.list(group_sizes), covariates = covs,
  design_formula = resolved$design_formula,
  organism = organism, gene_id_type = det$id_type,
  library_size_min = min(lib), library_size_median = median(lib), library_size_max = max(lib)
)
finish()
