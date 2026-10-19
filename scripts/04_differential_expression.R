#!/usr/bin/env Rscript
# =============================================================================
# 04_differential_expression.R  --  DESeq2 differential expression
# =============================================================================
# What it does
#   1. Builds the DESeq2 model  ~ covariates + condition  with the reference
#      group as the baseline.
#   2. Removes low-count genes (>= min_count reads in >= min_samples samples).
#   3. Runs DESeq() (size factors, dispersions, negative binomial GLM, Wald test).
#   4. Extracts test-vs-reference results with Benjamini-Hochberg adjusted
#      p-values (DESeq2 independent filtering on, alpha = padj_cutoff).
#   5. Shrinks log2 fold changes (apeglm by default) for ranking and plotting.
#      p-values and padj always come from the unshrunken Wald test.
#   6. Calls a gene significant when padj < padj_cutoff AND
#      |shrunken log2FC| > log2fc_cutoff.
#   7. Adds gene symbols and names (human / mouse) and per-group mean
#      normalised counts.
#
# Sign convention: positive log2FoldChange = higher in test_level than in
# reference_level.
#
# Outputs (in <run_dir>/04_differential_expression/)
#   de_results_all.csv          every gene tested, sorted by padj
#   de_results_significant.csv  significant genes only
#   normalized_counts.csv       DESeq2 size-factor normalised counts
#   de_summary.json             numbers, thresholds, design, method (read by the agent)
#   dispersion_plot.png         DESeq2 dispersion estimates (model-fit check)
#   dds.rds                     the fitted DESeqDataSet (for reruns / custom work)
#
# Usage
#   Rscript scripts/04_differential_expression.R --run_dir runs/my_run
#   Rscript scripts/04_differential_expression.R --run_dir runs/my_run --padj_cutoff 0.1 --log2fc_cutoff 0.58
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "scripts/x"), ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
})

STEP <- "04_differential_expression"
started <- Sys.time()

option_list <- list(
  make_option("--run_dir", type = "character", help = "Run folder (required)"),
  make_option("--config", type = "character", help = "Parameter file (default: <run_dir>/params.yaml)"),
  make_option("--set", type = "character", help = "Extra overrides"),
  make_option("--padj_cutoff", type = "double", help = "Adjusted p-value threshold (default 0.05)"),
  make_option("--log2fc_cutoff", type = "double", help = "|log2 fold change| threshold (default 1)"),
  make_option("--lfc_shrinkage", type = "character", help = "apeglm | ashr | normal | none"),
  make_option("--min_count", type = "integer", help = "Low-count filter: minimum reads"),
  make_option("--min_samples", type = "integer", help = "Low-count filter: minimum samples")
)
opt <- parse_args(OptionParser(option_list = option_list, description = "Step 04: DESeq2 differential expression."))
params <- load_params(opt)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
log_info("Step 04: differential expression with DESeq2 ", as.character(packageVersion("DESeq2")))
set.seed(params$seed)

inp <- read_validated_inputs(run_dir)
counts <- inp$counts
meta <- inp$meta
cond <- params$condition; ref <- params$reference_level; test <- params$test_level
covs <- params$covariates

# ---- Model data with R-safe names --------------------------------------------
# DESeq2 needs syntactic column and level names. We convert them here and keep
# a map so that every output uses the original names.
safe_name <- function(x) make.names(x)
col_data <- data.frame(row.names = meta$sample_id)
level_maps <- list()
for (v in c(covs, cond)) {
  sv <- safe_name(v)
  if (v %in% params$numeric_covariates) {
    col_data[[sv]] <- as.numeric(scale(as.numeric(meta[[v]])))   # centred and scaled
  } else {
    lv <- if (v == cond) levels(condition_factor(meta, params)) else sort(unique(as.character(meta[[v]])))
    map <- setNames(make.names(lv, unique = TRUE), lv)
    level_maps[[v]] <- as.list(map)
    col_data[[sv]] <- factor(unname(map[as.character(meta[[v]])]), levels = unname(map))
  }
}
cond_s <- safe_name(cond)
ref_s <- level_maps[[cond]][[ref]]; test_s <- level_maps[[cond]][[test]]
col_data[[cond_s]] <- relevel(col_data[[cond_s]], ref = ref_s)
design <- as.formula(paste("~", paste(c(safe_name(covs), cond_s), collapse = " + ")))
log_info("Design: ", params$design_formula, "   |   comparison: ", test, " vs ", ref)

# ---- Filter and fit ----------------------------------------------------------
keep <- low_count_filter(counts, params$min_count, params$min_samples)
log_info(sprintf("Low-count filter (>= %d reads in >= %d samples): %d of %d genes kept.",
                 params$min_count, params$min_samples, sum(keep), length(keep)))
dds <- DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE], colData = col_data, design = design)
dds <- with_logged_warnings(DESeq(dds, quiet = TRUE))
fit_type <- attr(dispersionFunction(dds), "fitType") %||% "parametric"
saveRDS(dds, file.path(out_dir, "dds.rds"))

coef_name <- paste(cond_s, test_s, "vs", ref_s, sep = "_")
if (!coef_name %in% resultsNames(dds)) {
  fail("Expected coefficient '", coef_name, "' not found among: ", paste(resultsNames(dds), collapse = ", "))
}
res <- results(dds, contrast = c(cond_s, test_s, ref_s), alpha = params$padj_cutoff)

# ---- Fold-change shrinkage ---------------------------------------------------
shrink_requested <- params$lfc_shrinkage %||% "apeglm"
shrink_used <- shrink_requested
lfc_shrunk <- res$log2FoldChange
lfc_se <- res$lfcSE
if (shrink_requested != "none") {
  shr <- tryCatch({
    if (shrink_requested == "apeglm") {
      if (!requireNamespace("apeglm", quietly = TRUE)) stop("package 'apeglm' is not installed")
      lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
    } else if (shrink_requested == "ashr") {
      if (!requireNamespace("ashr", quietly = TRUE)) stop("package 'ashr' is not installed")
      lfcShrink(dds, contrast = c(cond_s, test_s, ref_s), type = "ashr", quiet = TRUE)
    } else if (shrink_requested == "normal") {
      lfcShrink(dds, contrast = c(cond_s, test_s, ref_s), type = "normal", quiet = TRUE)
    } else stop("unknown lfc_shrinkage '", shrink_requested, "'")
  }, error = function(e) e)
  if (inherits(shr, "error")) {
    log_warn("Fold-change shrinkage '", shrink_requested, "' failed (", conditionMessage(shr),
             "). Falling back to 'normal' shrinkage.")
    shr <- lfcShrink(dds, contrast = c(cond_s, test_s, ref_s), type = "normal", quiet = TRUE)
    shrink_used <- "normal"
  }
  stopifnot(identical(rownames(shr), rownames(res)))
  lfc_shrunk <- shr$log2FoldChange
  lfc_se <- shr$lfcSE
}

# ---- Results table -----------------------------------------------------------
norm <- counts(dds, normalized = TRUE)
grp_raw <- as.character(meta[[cond]])
tab <- data.frame(
  gene_id = rownames(res),
  baseMean = res$baseMean,
  log2FoldChange = lfc_shrunk,
  lfcSE = lfc_se,
  log2FoldChange_unshrunken = res$log2FoldChange,
  stat = res$stat,
  pvalue = res$pvalue,
  padj = res$padj,
  stringsAsFactors = FALSE
)
tab[[paste0("mean_normalized_", ref)]] <- rowMeans(norm[, grp_raw == ref, drop = FALSE])
tab[[paste0("mean_normalized_", test)]] <- rowMeans(norm[, grp_raw == test, drop = FALSE])

sig <- !is.na(tab$padj) & tab$padj < params$padj_cutoff & abs(tab$log2FoldChange) > params$log2fc_cutoff
tab$significant <- sig
tab$direction <- ifelse(sig & tab$log2FoldChange > 0, "up",
                        ifelse(sig & tab$log2FoldChange < 0, "down", "not significant"))
tab$note <- ifelse(is.na(tab$pvalue), "no p-value: flagged as a count outlier (Cook's distance) by DESeq2",
                   ifelse(is.na(tab$padj), "no adjusted p-value: mean count too low (removed by independent filtering)", ""))

ann <- annotate_genes(tab$gene_id, params$organism, params$gene_id_type)
tab <- cbind(tab[, "gene_id", drop = FALSE], ann[, c("gene_symbol", "gene_name", "entrez_id")],
             tab[, setdiff(names(tab), "gene_id")])
tab <- tab[order(is.na(tab$padj), tab$padj, tab$pvalue), ]
num_cols <- c("baseMean", "log2FoldChange", "lfcSE", "log2FoldChange_unshrunken", "stat",
              paste0("mean_normalized_", c(ref, test)))
for (cc in num_cols) tab[[cc]] <- signif(tab[[cc]], 6)

write_table(tab, file.path(out_dir, "de_results_all.csv"))
write_table(tab[tab$significant, ], file.path(out_dir, "de_results_significant.csv"))

norm_out <- data.frame(gene_id = rownames(norm),
                       gene_symbol = ann$gene_symbol[match(rownames(norm), ann$gene_id)],
                       round(norm, 2), check.names = FALSE)
write_table(norm_out, file.path(out_dir, "normalized_counts.csv"))

png(file.path(out_dir, "dispersion_plot.png"), width = 7, height = 5.5, units = "in", res = 300)
plotDispEsts(dds, main = "DESeq2 dispersion estimates")
dev.off()
log_info("Saved plot: dispersion_plot.png")

# ---- Summary -----------------------------------------------------------------
n_up <- sum(tab$direction == "up"); n_down <- sum(tab$direction == "down")
top_genes <- function(dirn, n = 10) {
  x <- tab[tab$direction == dirn, ][seq_len(min(n, sum(tab$direction == dirn))), ]
  if (!nrow(x)) return(list())
  lapply(seq_len(nrow(x)), function(i) list(gene_id = x$gene_id[i], gene_symbol = x$gene_symbol[i],
                                            log2FoldChange = round(x$log2FoldChange[i], 3),
                                            padj = signif(x$padj[i], 3)))
}
summary <- list(
  method = sprintf("DESeq2 %s, Wald test", as.character(packageVersion("DESeq2"))),
  design_formula = params$design_formula,
  comparison = sprintf("%s vs %s", test, ref),
  sign_convention = sprintf("positive log2FoldChange = higher in '%s' than in '%s'", test, ref),
  condition = cond, reference_level = ref, test_level = test, covariates = covs,
  numeric_covariates_scaled = params$numeric_covariates,
  coefficient = coef_name,
  low_count_filter = sprintf(">= %d reads in >= %d samples", params$min_count, params$min_samples),
  n_genes_input = nrow(counts),
  n_genes_tested = nrow(tab),
  n_padj_na = sum(is.na(tab$padj)),
  n_cooks_outliers = sum(is.na(tab$pvalue)),
  padj_cutoff = params$padj_cutoff,
  log2fc_cutoff = params$log2fc_cutoff,
  significance_rule = sprintf("padj < %g and |shrunken log2FC| > %g", params$padj_cutoff, params$log2fc_cutoff),
  n_padj_below_cutoff = sum(tab$padj < params$padj_cutoff, na.rm = TRUE),
  n_significant = sum(tab$significant),
  n_up = n_up, n_down = n_down,
  lfc_shrinkage_requested = shrink_requested,
  lfc_shrinkage_used = shrink_used,
  p_adjust_method = "Benjamini-Hochberg (DESeq2 independent filtering on)",
  dispersion_fit_type = fit_type,
  size_factors = as.list(round(sizeFactors(dds), 4)),
  level_name_map = level_maps,
  top_up = top_genes("up"),
  top_down = top_genes("down")
)
write_json_file(summary, file.path(out_dir, "de_summary.json"))
log_info(sprintf("%d genes tested; %d significant (%d up, %d down in %s vs %s).",
                 nrow(tab), sum(tab$significant), n_up, n_down, test, ref))

write_step_provenance(run_dir, STEP, started,
                      params_used = params[c("condition", "reference_level", "test_level", "covariates",
                                             "design_formula", "min_count", "min_samples", "padj_cutoff",
                                             "log2fc_cutoff", "lfc_shrinkage", "organism", "seed")],
                      packages = c("DESeq2", "apeglm", "ashr", "AnnotationDbi", orgdb_for(params$organism)),
                      inputs = file.path(run_dir, "01_validate_inputs", c("counts_validated.csv", "metadata_validated.csv")),
                      outputs = list.files(out_dir),
                      notes = list(lfc_shrinkage_used = shrink_used))
log_info("Step 04 done.")
