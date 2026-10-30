#!/usr/bin/env Rscript
# =============================================================================
# I01_validate_isoforms.R  --  check the transcript quantification and design
# =============================================================================
# Checks the sample table, finds each sample's quantification (Salmon quant.sf
# or kallisto abundance.tsv), maps transcripts to genes, checks the design, and
# imports everything with tximport (transcript level):
#   txi_counts.rds      countsFromAbundance = "no"   -> DEI (DESeq2 with length offsets)
#   txi_dtu_scaled.rds  countsFromAbundance = "dtuScaledTPM" -> DUI (DEXSeq)
#
# Outputs (in <run_dir>/I01_validate_isoforms/): validation_report.json / .md,
# samples_validated.csv, transcripts.csv, txi_*.rds, resolved_params.yaml
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

STEP <- "I01_validate_isoforms"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "Isoform step I01: check the inputs."))
if (is_blank(opt$run_dir)) fail("--run_dir is required.")
params <- load_isoform_params(opt$run_dir, use_resolved = FALSE)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
log_info("Step I01: checking the transcript quantification and the design")

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
  writeLines(c("# Isoform input check", "", sprintf("Overall: **%s**", overall), "",
               vapply(checks, function(c) sprintf("- **%s** (%s): %s", c$id, c$status, c$message), "")),
             file.path(out_dir, "validation_report.md"))
  overall
}
stop_now <- function(msg) { finish(); fail(msg) }
for (pkg in c("tximport", "DESeq2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) { check("packages", "fail", sprintf("Package '%s' is not installed.", pkg)); stop_now("Missing package.") }
}

# ---- Sample table -----------------------------------------------------------------
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
n_ref <- sum(meta[[cond]] == ref); n_test <- sum(meta[[cond]] == test)
if (min(n_ref, n_test) < 2) { check("replicates", "fail", sprintf("%s has %d and %s has %d sample(s); both need at least 2.", ref, n_ref, test, n_test)); stop_now("Too few replicates.") }
check("replicates", if (min(n_ref, n_test) < 3) "warning" else "pass", sprintf("%s: %d samples, %s: %d samples%s.", ref, n_ref, test, n_test,
                                                                                 if (min(n_ref, n_test) < 3) "; at least 3 per group are recommended" else ""))
covs <- params$covariates
miss <- setdiff(covs, names(meta))
if (length(miss)) { check("covariates", "fail", paste("Covariate column(s) not found:", paste(miss, collapse = ", "))); stop_now("Missing covariates.") }
design_ok <- function(m) {
  if (!nrow(m)) return(FALSE)
  df <- as.data.frame(lapply(m[, c(covs, cond), drop = FALSE], factor))
  X <- model.matrix(as.formula(paste("~", paste(make.names(c(covs, cond)), collapse = " + "))), setNames(df, make.names(names(df))))
  qr(X)$rank == ncol(X) && nrow(X) > ncol(X)
}
if (!design_ok(meta)) { check("design_full_rank", "fail", sprintf("The design ~ %s cannot be estimated: a covariate lines up with '%s', or there are too few samples.", paste(c(covs, cond), collapse = " + "), cond)); stop_now("Design not estimable.") }
pair <- meta[meta[[cond]] %in% c(ref, test), ]
if (!design_ok(pair)) { check("design_full_rank", "fail", "Within the two compared groups the covariates cannot be separated from the comparison (needed for isoform usage)."); stop_now("Design not estimable.") }
check("design_full_rank", "pass", sprintf("Design ~ %s can be estimated.", paste(c(covs, cond), collapse = " + ")))

# ---- Quantification files --------------------------------------------------------------
qdir <- resolve_path(params$quant_dir)
if (is.null(qdir) || !dir.exists(qdir)) { check("quant_dir", "fail", "quant_dir was not found."); stop_now("No quantification folder.") }
qfile <- quant_file_name(params$quant_type)
files <- setNames(file.path(qdir, meta$sample_id, qfile), meta$sample_id)
if (!all(file.exists(files))) {
  check("quant_files", "fail", sprintf("No %s for sample(s): %s (expected <quant_dir>/<sample>/%s).", qfile,
                                       paste(names(files)[!file.exists(files)], collapse = ", "), qfile))
  stop_now("Missing quantification files.")
}
check("quant_files", "pass", sprintf("%s found for all %d samples.", qfile, length(files)))
# copy into the run so it is self-contained
dir.create(file.path(run_dir, "inputs", "quant"), recursive = TRUE, showWarnings = FALSE)
for (s in names(files)) {
  dir.create(file.path(run_dir, "inputs", "quant", s), showWarnings = FALSE)
  file.copy(files[[s]], file.path(run_dir, "inputs", "quant", s, qfile), overwrite = TRUE)
}
files <- setNames(file.path(run_dir, "inputs", "quant", names(files), qfile), names(files))

# ---- Transcript -> gene ----------------------------------------------------------------
tx2gene <- if (!is_blank(params$tx2gene_file)) {
  t <- read_table_file(resolve_path(params$tx2gene_file))
  if (!all(c("transcript_id", "gene_id") %in% names(t))) { check("tx2gene", "fail", "tx2gene_file needs columns transcript_id and gene_id."); stop_now("Bad tx2gene.") }
  if (!"gene_name" %in% names(t)) t$gene_name <- NA_character_
  t[, c("transcript_id", "gene_id", "gene_name")]
} else if (!is_blank(params$gtf_file) && file.exists(resolve_path(params$gtf_file))) {
  tx2gene_from_gtf(resolve_path(params$gtf_file))
} else { check("tx2gene", "fail", "Give gtf_file or tx2gene_file (transcript -> gene)."); stop_now("No annotation.") }
if (isTRUE(params$ignore_tx_version)) tx2gene$transcript_id <- sub("\\.\\d+$", "", tx2gene$transcript_id)
first <- data.table::fread(files[[1]], data.table = FALSE, select = 1)[[1]]
strip_bar <- grepl("\\|", first[1])
names_q <- if (strip_bar) sub("\\|.*$", "", first) else first
if (isTRUE(params$ignore_tx_version)) names_q <- sub("\\.\\d+$", "", names_q)
matched <- mean(names_q %in% tx2gene$transcript_id)
if (matched < 0.9) {
  check("tx2gene", if (matched < 0.5) "fail" else "warning",
        sprintf("Only %.0f%% of the quantified transcripts are in the annotation. Use the GTF the index was built from, or ignore_tx_version: true.", 100 * matched))
  if (matched < 0.5) stop_now("Annotation does not match.")
} else check("tx2gene", "pass", sprintf("%.1f%% of %d transcripts mapped to %d genes.", 100 * matched, length(names_q),
                                        length(unique(tx2gene$gene_id[tx2gene$transcript_id %in% names_q]))))

# ---- Import ------------------------------------------------------------------------------
imp <- function(cfa) suppressMessages(tximport::tximport(files, type = params$quant_type, txOut = TRUE, countsFromAbundance = cfa,
                                                        tx2gene = tx2gene[, c("transcript_id", "gene_id")],
                                                        ignoreTxVersion = isTRUE(params$ignore_tx_version), ignoreAfterBar = strip_bar))
txi <- imp("no")
txi_dtu <- imp("dtuScaledTPM")
clean <- function(x) { if (strip_bar) x <- sub("\\|.*$", "", x); if (isTRUE(params$ignore_tx_version)) x <- sub("\\.\\d+$", "", x); x }
for (k in c("abundance", "counts", "length")) { rownames(txi[[k]]) <- clean(rownames(txi[[k]])); rownames(txi_dtu[[k]]) <- clean(rownames(txi_dtu[[k]])) }
keep_tx <- rownames(txi$counts) %in% tx2gene$transcript_id
for (k in c("abundance", "counts", "length")) { txi[[k]] <- txi[[k]][keep_tx, , drop = FALSE]; txi_dtu[[k]] <- txi_dtu[[k]][keep_tx, , drop = FALSE] }
saveRDS(txi, file.path(out_dir, "txi_counts.rds"))
saveRDS(txi_dtu, file.path(out_dir, "txi_dtu_scaled.rds"))
tx <- tx2gene[match(rownames(txi$counts), tx2gene$transcript_id), ]
tx$n_transcripts_in_gene <- as.integer(table(tx$gene_id)[tx$gene_id])
write_table(tx, file.path(out_dir, "transcripts.csv"))
write_table(meta, file.path(out_dir, "samples_validated.csv"))
lib <- colSums(txi$counts)
low <- names(lib)[lib < 1e6]
check("library_sizes", if (length(low)) "warning" else "pass",
      if (length(low)) sprintf("Fewer than 1 million assigned reads: %s.", paste(low, collapse = ", ")) else
        sprintf("Assigned reads per sample: %s to %s.", format(round(min(lib)), big.mark = ","), format(round(max(lib)), big.mark = ",")))
organism <- params$organism %||% "auto"
if (organism == "auto") organism <- if (mean(grepl("^ENSMUSG", head(tx$gene_id, 500))) > 0.5) "mouse" else if (mean(grepl("^ENSG\\d", head(tx$gene_id, 500))) > 0.5) "human" else "other"

if (finish() == "fail") fail("The input check failed.")
write_yaml_file(list(reference_level = ref, test_level = test, organism = organism,
                     min_samples = params$min_samples %||% min(n_ref, n_test), n_smallest_group = min(n_ref, n_test),
                     design_formula = paste("~", paste(c(covs, cond), collapse = " + ")),
                     n_samples = nrow(meta), n_transcripts = nrow(txi$counts), n_genes = length(unique(tx$gene_id)),
                     group_sizes = setNames(as.list(as.integer(table(meta[[cond]]))), names(table(meta[[cond]])))),
                file.path(out_dir, "resolved_params.yaml"), header = "Values used by all later isoform steps (filled in by I01).")
write_step_provenance(run_dir, STEP, started, params_used = params[c("quant_type", "condition", "covariates", "ignore_tx_version")],
                      packages = c("tximport"), inputs = unname(files), outputs = list.files(out_dir))
log_info(sprintf("Step I01 done: %d transcripts in %d genes, %d samples.", nrow(txi$counts), length(unique(tx$gene_id)), nrow(meta)))
