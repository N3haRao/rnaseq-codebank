#!/usr/bin/env Rscript
# =============================================================================
# I02_differential_isoform_expression.R  --  DEI: DESeq2 per transcript
# =============================================================================
# The main pipeline's differential-expression method applied to transcripts:
# tximport transcript counts with effective-length offsets
# (DESeqDataSetFromTximport), the same low-count filter, design (covariates +
# condition), Wald test, BH adjustment and fold-change shrinkage as step 04.
# A transcript is significant when padj < padj_cutoff and |shrunken log2FC| >
# log2fc_cutoff. Positive log2FC = higher in the test group.
#
# Outputs (in <run_dir>/I02_differential_isoform_expression/):
#   dei_results_all.csv, dei_results_significant.csv, dei_summary.json,
#   dei_volcano.png, dei_ma_plot.png, dds_transcripts.rds
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "isoforms/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "R", "gtf.R"))
source(file.path(REPO_ROOT, "isoforms", "R", "isoform_utils.R"))
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(DESeq2)
})

STEP <- "I02_differential_isoform_expression"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "Isoform step I02: differential isoform expression."))
if (is_blank(opt$run_dir)) fail("--run_dir is required.")
params <- load_isoform_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
set.seed(params$seed)
log_info("Step I02: differential isoform expression with DESeq2 ", as.character(packageVersion("DESeq2")))

inp <- read_isoform_inputs(run_dir)
meta <- inp$meta; txi <- inp$txi; tx <- inp$tx
cond <- params$condition; ref <- params$reference_level; test <- params$test_level
mf <- model_frame(meta, params)
cond_s <- make.names(cond)
ref_s <- mf$maps[[cond]][[ref]]; test_s <- mf$maps[[cond]][[test]]
design <- as.formula(paste("~", paste(make.names(c(params$covariates, cond)), collapse = " + ")))

keep <- low_count_filter(round(txi$counts), params$min_count, params$min_samples)
log_info(sprintf("Low-count filter (>= %d reads in >= %d samples): %d of %d transcripts kept.",
                 params$min_count, params$min_samples, sum(keep), length(keep)))
for (k in c("abundance", "counts", "length")) txi[[k]] <- txi[[k]][keep, , drop = FALSE]
dds <- suppressMessages(DESeqDataSetFromTximport(txi, colData = mf$col_data, design = design))
dds <- with_logged_warnings(DESeq(dds, quiet = TRUE))
saveRDS(dds, file.path(out_dir, "dds_transcripts.rds"))
coef_name <- paste(cond_s, test_s, "vs", ref_s, sep = "_")
res <- results(dds, contrast = c(cond_s, test_s, ref_s), alpha = params$padj_cutoff)

shrink_used <- params$lfc_shrinkage %||% "apeglm"
lfc <- res$log2FoldChange; lfc_se <- res$lfcSE
if (shrink_used != "none") {
  shr <- tryCatch(switch(shrink_used,
    apeglm = lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE),
    ashr = lfcShrink(dds, contrast = c(cond_s, test_s, ref_s), type = "ashr", quiet = TRUE),
    normal = lfcShrink(dds, contrast = c(cond_s, test_s, ref_s), type = "normal", quiet = TRUE),
    stop("unknown lfc_shrinkage")), error = function(e) e)
  if (inherits(shr, "error")) {
    log_warn("Shrinkage '", shrink_used, "' failed (", conditionMessage(shr), "); using 'normal'.")
    shr <- lfcShrink(dds, contrast = c(cond_s, test_s, ref_s), type = "normal", quiet = TRUE)
    shrink_used <- "normal"
  }
  lfc <- shr$log2FoldChange; lfc_se <- shr$lfcSE
}

grp <- meta[[cond]]
tpm <- txi$abundance
tab <- data.frame(transcript_id = rownames(res), tx[match(rownames(res), tx$transcript_id), c("gene_id", "gene_name", "n_transcripts_in_gene")],
                  baseMean = signif(res$baseMean, 6), log2FoldChange = signif(lfc, 6), lfcSE = signif(lfc_se, 6),
                  log2FoldChange_unshrunken = signif(res$log2FoldChange, 6), stat = res$stat, pvalue = res$pvalue, padj = res$padj,
                  check.names = FALSE, stringsAsFactors = FALSE)
tab[[paste0("mean_TPM_", ref)]] <- signif(rowMeans(tpm[, grp == ref, drop = FALSE]), 5)
tab[[paste0("mean_TPM_", test)]] <- signif(rowMeans(tpm[, grp == test, drop = FALSE]), 5)
tab$significant <- !is.na(tab$padj) & tab$padj < params$padj_cutoff & abs(tab$log2FoldChange) > params$log2fc_cutoff
tab$direction <- ifelse(tab$significant & tab$log2FoldChange > 0, "up", ifelse(tab$significant & tab$log2FoldChange < 0, "down", "not significant"))
tab$note <- ifelse(is.na(tab$pvalue), "no p-value: count outlier (Cook's distance)",
                   ifelse(is.na(tab$padj), "no adjusted p-value: removed by independent filtering", ""))
tab <- tab[order(is.na(tab$padj), tab$padj, tab$pvalue), ]
write_table(tab, file.path(out_dir, "dei_results_all.csv"))
sig <- tab[tab$significant, ]
write_table(sig, file.path(out_dir, "dei_results_significant.csv"))
rule <- sprintf("padj < %g and |shrunken log2FC| > %g", params$padj_cutoff, params$log2fc_cutoff)
log_info(sprintf("%d significant transcripts (%d up, %d down) in %d genes; %s.", nrow(sig), sum(sig$direction == "up"),
                 sum(sig$direction == "down"), length(unique(sig$gene_id)), rule))

# ---- Figures -----------------------------------------------------------------
dir_cols <- direction_colours(params)
lab_of <- function(d) ifelse(is.na(d$gene_name) | d$gene_name == "", d$transcript_id, paste0(d$gene_name, " (", d$transcript_id, ")"))
v <- tab[!is.na(tab$pvalue), ]
v$neglog10p <- -log10(pmax(v$pvalue, min(v$pvalue[v$pvalue > 0], 1e-300)))
v$direction <- factor(v$direction, levels = names(dir_cols))
v <- v[order(v$direction != "not significant"), ]
labs_df <- head(v[v$significant, ][order(v$padj[v$significant]), ], params$n_label)
labs_df$lab <- lab_of(labs_df)
p <- ggplot(v, aes(log2FoldChange, neglog10p, colour = direction)) + geom_point(size = 0.9, alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1) * params$log2fc_cutoff, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = dir_cols, drop = FALSE, name = NULL,
                      labels = c(up = sprintf("Higher in %s (%d)", test, sum(sig$direction == "up")),
                                 down = sprintf("Lower in %s (%d)", test, sum(sig$direction == "down")), `not significant` = "Not significant")) +
  labs(title = sprintf("Differential isoform expression: %s vs %s", test, ref), subtitle = paste("Significant:", rule),
       x = sprintf("log2 fold change (%s / %s, shrunken)", test, ref), y = "-log10(p-value)") +
  theme_codebank() + theme(legend.position = "bottom")
if (nrow(labs_df) && requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(data = labs_df, aes(label = lab), colour = "black", size = 2.6, max.overlaps = 40, show.legend = FALSE)
}
save_plot(p, file.path(out_dir, "dei_volcano.png"), width = 8, height = 6.5)
m <- tab[tab$baseMean > 0, ]
m$direction <- factor(m$direction, levels = names(dir_cols))
p_ma <- ggplot(m[order(m$direction != "not significant"), ], aes(log10(baseMean), log2FoldChange, colour = direction)) +
  geom_point(size = 0.8, alpha = 0.7) + geom_hline(yintercept = 0, colour = "grey30") +
  scale_colour_manual(values = dir_cols, drop = FALSE, name = NULL) +
  labs(title = sprintf("MA plot, transcripts: %s vs %s", test, ref), x = "log10 mean normalised count", y = "log2 fold change (shrunken)") +
  theme_codebank() + theme(legend.position = "bottom")
save_plot(p_ma, file.path(out_dir, "dei_ma_plot.png"), width = 7.5, height = 5.5)

top <- function(dirn, n = 10) {
  x <- head(sig[sig$direction == dirn, ], n)
  if (!nrow(x)) return(list())
  lapply(seq_len(nrow(x)), function(i) list(transcript_id = x$transcript_id[i], gene_id = x$gene_id[i], gene_name = x$gene_name[i],
                                            log2FoldChange = round(x$log2FoldChange[i], 3), padj = signif(x$padj[i], 3)))
}
write_json_file(list(
  method = sprintf("DESeq2 %s on transcript counts (tximport, effective-length offsets), Wald test", as.character(packageVersion("DESeq2"))),
  design_formula = params$design_formula, comparison = sprintf("%s vs %s", test, ref),
  sign_convention = sprintf("positive log2FoldChange = higher in '%s'", test),
  low_count_filter = sprintf(">= %d reads in >= %d samples", params$min_count, params$min_samples),
  significance_rule = rule, lfc_shrinkage_used = shrink_used,
  n_transcripts_tested = sum(!is.na(tab$padj)), n_significant = nrow(sig), n_up = sum(sig$direction == "up"),
  n_down = sum(sig$direction == "down"), n_genes_with_significant_transcripts = length(unique(sig$gene_id)),
  top_up = top("up"), top_down = top("down")
), file.path(out_dir, "dei_summary.json"))
write_step_provenance(run_dir, STEP, started,
                      params_used = params[c("min_count", "min_samples", "lfc_shrinkage", "padj_cutoff", "log2fc_cutoff", "seed")],
                      packages = c("DESeq2", "apeglm", "tximport", "ggplot2"), outputs = list.files(out_dir))
log_info("Step I02 done.")
