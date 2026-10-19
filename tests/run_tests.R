#!/usr/bin/env Rscript
# =============================================================================
# tests/run_tests.R  --  end-to-end tests of the code bank
# =============================================================================
# Run from the repository root:
#   Rscript tests/run_tests.R                 # everything
#   Rscript tests/run_tests.R --skip_enrichment   # if offline / no MSigDB download
#
# What is tested
#   A. airway: the full pipeline runs, every expected file exists
#   B. airway: DE results are IDENTICAL to a hand-written DESeq2 analysis with
#      the same settings (log2FC, p-values, padj)
#   C. airway: well-known dexamethasone-induced genes are significant and up
#      (DUSP1, KLF15, PER1, TSC22D3, FKBP5, ZBTB16)
#   D. airway: a rerun from the saved parameters gives identical results; exclude_samples works
#   E. pasilla: covariate design, non-human organism, enrichment skipped cleanly
#   F. planted errors: step 01 fails with the right check for each broken input
#   G. swapped sample labels: step 03 flags the swapped samples
#
# Test runs are written to tests/output/ (ignored by git). Exit code 0 = all passed.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "tests/x"), ".."))
})
setwd(REPO_ROOT)
source(file.path(REPO_ROOT, "R", "utils.R"))
args <- commandArgs(trailingOnly = TRUE)
skip_enrichment <- "--skip_enrichment" %in% args

out_root <- file.path(REPO_ROOT, "tests", "output")
unlink(out_root, recursive = TRUE)
dir.create(out_root, recursive = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")

results <- data.frame(test = character(0), status = character(0), detail = character(0))
record <- function(test, ok, detail = "") {
  results[nrow(results) + 1, ] <<- list(test, if (isTRUE(ok)) "PASS" else "FAIL", detail)
  cat(sprintf("  [%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", test, if (nzchar(detail)) paste0("- ", detail) else ""))
}

# Run the pipeline; returns list(code, run_dir)
run_pipe <- function(args_vec, runs_dir = out_root) {
  before <- list.dirs(runs_dir, recursive = FALSE)
  log <- tempfile(fileext = ".log")
  code <- system2(rscript, c(file.path("scripts", "run_pipeline.R"), "--runs_dir", shQuote(runs_dir), args_vec),
                  stdout = log, stderr = log)
  after <- setdiff(list.dirs(runs_dir, recursive = FALSE), before)
  list(code = code, run_dir = if (length(after)) after[1] else NA, log = log)
}

enrich_set <- if (skip_enrichment) "run_enrichment=false" else NULL
set_arg <- function(...) { s <- c(...); if (length(s)) c("--set", shQuote(paste(s, collapse = ";"))) else character(0) }

# ---- Manifest consistency -----------------------------------------------------
cat("\nManifest\n")
man <- yaml::read_yaml("manifest/tools.yaml")
defs <- load_defaults()
man_params <- unique(unlist(lapply(man$steps, function(s) names(s$parameters))))
record("every manifest parameter has a default", all(man_params %in% names(defs)),
       paste(setdiff(man_params, names(defs)), collapse = ", "))
record("every manifest step has a script", all(file.exists(vapply(man$steps, `[[`, "", "script"))))
record("manifest steps match the pipeline order",
       identical(vapply(man$steps, `[[`, "", "name"), vapply(PIPELINE_STEPS, `[[`, "", "name")))

# ---- Example data -------------------------------------------------------------
cat("\nPreparing example data\n")
need <- c("data/example/airway/airway_counts.csv", "data/example/pasilla/pasilla_counts.tsv")
if (!all(file.exists(need))) system2(rscript, "scripts/make_example_data.R")
record("example data present", all(file.exists(need)))

# ---- A. airway full run ---------------------------------------------------------
cat("\nA. airway full pipeline\n")
a <- run_pipe(c("--config", "config/examples/airway.yaml", set_arg("run_name=test_airway", enrich_set)))
record("pipeline exit code 0", a$code == 0, a$run_dir)
rd <- a$run_dir
expected <- c("01_validate_inputs/validation_report.json", "01_validate_inputs/resolved_params.yaml",
              "02_qc/qc_metrics.csv", "02_qc/library_sizes.png",
              "03_exploratory/pca_plot.png", "03_exploratory/outlier_flags.csv", "03_exploratory/vsd.rds",
              "04_differential_expression/de_results_all.csv", "04_differential_expression/de_summary.json",
              "05_de_plots/volcano_plot.png", "05_de_plots/ma_plot.png", "05_de_plots/heatmap_top_degs.png",
              "06_enrichment/enrichment_summary.json", "07_report/methods.md",
              "provenance.json", "run_summary.json")
missing <- expected[!file.exists(file.path(rd, expected))]
record("expected output files exist", !length(missing), paste(missing, collapse = ", "))
summ <- read_json_file(file.path(rd, "run_summary.json"))
if (!skip_enrichment) {
  en <- read_json_file(file.path(rd, "06_enrichment", "enrichment_summary.json"))
  record("enrichment ran (needs internet on first use)", identical(en$status, "done"),
         if (!identical(en$status, "done")) "see logs/06_enrichment.log" else "")
}
if (rmarkdown::pandoc_available()) record("HTML report written", file.exists(file.path(rd, "07_report/report.html")))

# ---- B. identical to a hand-written DESeq2 analysis -----------------------------
cat("\nB. comparison with a hand-written DESeq2 analysis\n")
suppressPackageStartupMessages(library(DESeq2))
cnt <- as.matrix(read.csv("data/example/airway/airway_counts.csv", row.names = 1, check.names = FALSE))
md <- read.csv("data/example/airway/airway_metadata.csv", row.names = 1)
md$condition <- factor(md$condition, levels = c("control", "dexamethasone"))
md$cell_line <- factor(md$cell_line)
keep <- rowSums(cnt >= 10) >= 4
dds <- DESeqDataSetFromMatrix(cnt[keep, ], md, ~ cell_line + condition)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", "dexamethasone", "control"), alpha = 0.05)
shr <- lfcShrink(dds, coef = "condition_dexamethasone_vs_control", type = "apeglm", quiet = TRUE)
mine <- read.csv(file.path(rd, "04_differential_expression", "de_results_all.csv"))
mine <- mine[match(rownames(res), mine$gene_id), ]
same <- function(x, y, tol) isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol))
record("same genes tested", identical(sort(mine$gene_id), sort(rownames(res))), sprintf("%d genes", nrow(res)))
record("p-values identical", same(mine$pvalue, res$pvalue, 1e-8))
record("adjusted p-values identical", same(mine$padj, res$padj, 1e-8))
record("unshrunken log2FC identical", same(mine$log2FoldChange_unshrunken, res$log2FoldChange, 1e-5))
record("apeglm log2FC identical", same(mine$log2FoldChange, shr$log2FoldChange, 1e-5))
n_sig <- sum(res$padj < 0.05 & abs(shr$log2FoldChange) > 1, na.rm = TRUE)
record("significant gene count identical", n_sig == summ$n_significant, sprintf("%d vs %s", n_sig, summ$n_significant))

# ---- C. known biology ----------------------------------------------------------------
cat("\nC. known dexamethasone-response genes\n")
known <- c(DUSP1 = "ENSG00000120129", KLF15 = "ENSG00000163884", PER1 = "ENSG00000179094",
           TSC22D3 = "ENSG00000157514", FKBP5 = "ENSG00000096060", ZBTB16 = "ENSG00000109906")
de_all <- read.csv(file.path(rd, "04_differential_expression", "de_results_all.csv"))
hit <- de_all[de_all$gene_id %in% known, c("gene_symbol", "direction", "padj")]
n_up <- sum(hit$direction == "up")
record("known glucocorticoid targets are significant and up", n_up >= 5,
       sprintf("%d of %d: %s", n_up, length(known), paste(hit$gene_symbol[hit$direction == "up"], collapse = ", ")))

# ---- D. reproducibility ---------------------------------------------------------------
cat("\nD. rerun from saved parameters\n")
d <- run_pipe(c("--based_on", shQuote(rd), "--steps", "1-4"))
r2 <- read.csv(file.path(d$run_dir, "04_differential_expression", "de_results_all.csv"))
r1 <- read.csv(file.path(rd, "04_differential_expression", "de_results_all.csv"))
record("rerun gives identical DE table", d$code == 0 && isTRUE(all.equal(r1, r2)))

x <- run_pipe(c("--based_on", shQuote(rd), "--steps", "1-4", set_arg("exclude_samples=[SRR1039508]")))
xs <- read_json_file(file.path(x$run_dir, "01_validate_inputs", "validation_report.json"))
record("exclude_samples leaves a sample out", x$code == 0 && xs$summary$n_samples == 7)

# ---- E. pasilla --------------------------------------------------------------------------
cat("\nE. pasilla (covariate, non-model organism)\n")
e <- run_pipe(c("--config", "config/examples/pasilla.yaml", set_arg("run_name=test_pasilla")))
record("pipeline exit code 0", e$code == 0, e$run_dir)
pde <- read_json_file(file.path(e$run_dir, "04_differential_expression", "de_summary.json"))
record("design includes the type covariate", identical(pde$design_formula, "~ type + condition"))
record("significant genes found", isTRUE(pde$n_significant > 50), sprintf("%s genes", pde$n_significant))
pen <- read_json_file(file.path(e$run_dir, "06_enrichment", "enrichment_summary.json"))
record("enrichment skipped with a reason", identical(pen$status, "skipped"))

# ---- F. planted errors ---------------------------------------------------------------------
cat("\nF. planted errors are caught by step 01\n")
bad_dir <- file.path(out_root, "bad_inputs"); dir.create(bad_dir)
cnt_df <- read.csv("data/example/airway/airway_counts.csv", check.names = FALSE)
meta_df <- read.csv("data/example/airway/airway_metadata.csv")
write_case <- function(name, counts = cnt_df, meta = meta_df) {
  cf <- file.path(bad_dir, paste0(name, "_counts.csv")); mf <- file.path(bad_dir, paste0(name, "_metadata.csv"))
  write.csv(counts, cf, row.names = FALSE); write.csv(meta, mf, row.names = FALSE)
  c(cf, mf)
}
expect_fail <- function(name, check_id, files, extra_set = NULL) {
  r <- run_pipe(c("--counts_file", shQuote(files[1]), "--metadata_file", shQuote(files[2]),
                  "--condition", "condition", "--reference_level", "control",
                  "--steps", "1", set_arg(paste0("run_name=bad_", name), extra_set)))
  rep <- read_json_file(file.path(r$run_dir, "01_validate_inputs", "validation_report.json"))
  failed <- rep$checks$id[rep$checks$status == "fail"]
  record(sprintf("%-28s -> fails at '%s'", name, check_id), r$code != 0 && check_id %in% failed,
         if (length(failed)) rep$checks$message[rep$checks$status == "fail"][1] else "no failure raised")
}
m1 <- meta_df; m1$sample_id[3] <- paste0(m1$sample_id[3], "_x")
expect_fail("mismatched sample name", "sample_matching", write_case("mismatch", meta = m1))
tpm <- cnt_df; tpm[-1] <- lapply(tpm[-1], function(x) x / sum(x) * 1e6)
expect_fail("TPM instead of counts", "count_values", write_case("tpm", counts = tpm))
lg <- cnt_df; lg[-1] <- lapply(lg[-1], function(x) log2(x / sum(x) * 1e6 + 1))
expect_fail("log-transformed values", "count_values", write_case("log", counts = lg))
na <- cnt_df; na[5, 3] <- NA
expect_fail("missing value", "missing_values", write_case("na", counts = na))
neg <- cnt_df; neg[5, 3] <- -4
expect_fail("negative value", "negative_values", write_case("negative", counts = neg))
m2 <- meta_df; names(m2)[names(m2) == "condition"] <- "treatment"
expect_fail("missing group column", "condition_column", write_case("nocol", meta = m2))
m3 <- meta_df; m3$batch <- m3$condition
expect_fail("batch confounded with group", "design_full_rank", write_case("confounded", meta = m3), "covariates=[batch]")
# 2 controls + 1 dexamethasone sample
expect_fail("one replicate in a group", "replicates", write_case("onerep", counts = cnt_df[, 1:4], meta = meta_df[1:3, ]))
expect_fail("wrong reference name", "reference_level", write_case("badref"), "reference_level=Control")

# ---- G. swapped labels -----------------------------------------------------------------------
cat("\nG. swapped labels are flagged in the PCA check\n")
m4 <- meta_df; i <- c(1, 2); m4$condition[i] <- rev(m4$condition[i])   # swap one control / treated pair
f4 <- write_case("swapped", meta = m4)
g <- run_pipe(c("--counts_file", shQuote(f4[1]), "--metadata_file", shQuote(f4[2]), "--condition", "condition",
                "--reference_level", "control", "--test_level", "dexamethasone", "--steps", "1-3",
                set_arg("run_name=swapped")))
fl <- read.csv(file.path(g$run_dir, "03_exploratory", "outlier_flags.csv"))
flagged <- fl$sample_id[fl$flag == "clusters_with_other_group"]
record("both swapped samples flagged", all(meta_df$sample_id[i] %in% flagged),
       paste("flagged:", paste(flagged, collapse = ", ")))

# ---- Result ------------------------------------------------------------------------------------
cat("\n================================================================\n")
n_fail <- sum(results$status == "FAIL")
cat(sprintf(" %d tests: %d passed, %d failed\n", nrow(results), sum(results$status == "PASS"), n_fail))
if (n_fail) print(results[results$status == "FAIL", ], row.names = FALSE)
cat("================================================================\n")
write.csv(results, file.path(out_root, "test_results.csv"), row.names = FALSE)
quit(save = "no", status = if (n_fail) 1 else 0)
