#!/usr/bin/env Rscript
# =============================================================================
# I03_differential_isoform_usage.R  --  DUI: DEXSeq on transcripts + stageR
# =============================================================================
# Does the SHARE of a gene's expression that each transcript makes up change
# between the groups (an isoform switch), regardless of the gene's total level?
# Method (Love, Soneson & Patro 2018, F1000Research "Swimming downstream"):
#   * transcript counts scaled for DTU (tximport countsFromAbundance = "dtuScaledTPM")
#   * filter as DRIMSeq::dmFilter: a transcript needs >= min_feature_count reads and
#     >= min_feature_prop of its gene in >= n samples of the smaller group; a gene
#     needs >= min_gene_count reads in every sample and >= 2 kept transcripts
#   * DEXSeq: ~ sample + exon (+ covariate:exon) + condition:exon vs the model
#     without condition:exon, transcripts as "exons" of their gene
#   * stageR two-stage testing: genes screened with DEXSeq's per-gene q-value,
#     then transcripts confirmed within significant genes; overall FDR = dui_alpha
# Only the two compared groups are used. A transcript is called when its gene
# passes the screen, the transcript is confirmed, and |delta usage| >= min_delta_usage.
# "up" = larger share of the gene in the test group.
#
# Outputs (in <run_dir>/I03_differential_isoform_usage/):
#   dui_transcripts_all.csv, dui_transcripts_significant.csv, dui_genes.csv,
#   dui_summary.json, dui_volcano.png, dui_top_genes_usage.png
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
})

STEP <- "I03_differential_isoform_usage"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "Isoform step I03: differential isoform usage."))
if (is_blank(opt$run_dir)) fail("--run_dir is required.")
params <- load_isoform_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
set.seed(params$seed)
if (!requireNamespace("DEXSeq", quietly = TRUE)) fail("Package 'DEXSeq' is not installed (BiocManager::install('DEXSeq')).")
use_stager <- requireNamespace("stageR", quietly = TRUE)
log_info("Step I03: differential isoform usage with DEXSeq ", as.character(packageVersion("DEXSeq")),
         if (use_stager) paste0(" and stageR ", packageVersion("stageR")) else " (stageR not installed: per-transcript BH instead)")

inp <- read_isoform_inputs(run_dir)
cond <- params$condition; ref <- params$reference_level; test <- params$test_level
meta <- inp$meta[inp$meta[[cond]] %in% c(ref, test), ]
cnt <- inp$txi_dtu$counts[, meta$sample_id, drop = FALSE]
tx <- inp$tx
gene_of <- tx$gene_id[match(rownames(cnt), tx$transcript_id)]
n_small <- min(table(meta[[cond]]))

# ---- Filter (DRIMSeq::dmFilter rules) ------------------------------------------------
share <- function(m, g) { tot <- rowsum(m, g, reorder = FALSE)[g, , drop = FALSE]; p <- m / tot; p[!is.finite(p)] <- 0; p }
prop <- share(cnt, gene_of)
keep <- rowSums(cnt >= params$min_feature_count) >= n_small & rowSums(prop >= params$min_feature_prop) >= n_small
cnt <- cnt[keep, , drop = FALSE]; gene_of <- gene_of[keep]
gene_tot <- rowsum(cnt, gene_of, reorder = FALSE)
n_tx <- table(gene_of)
ok_genes <- rownames(gene_tot)[apply(gene_tot >= params$min_gene_count, 1, all) & as.integer(n_tx[rownames(gene_tot)]) >= 2]
sel <- gene_of %in% ok_genes
cnt <- cnt[sel, , drop = FALSE]; gene_of <- gene_of[sel]
log_info(sprintf("After filtering: %d transcripts in %d genes with >= 2 expressed transcripts (from %d transcripts).",
                 nrow(cnt), length(unique(gene_of)), nrow(inp$txi_dtu$counts)))
if (length(unique(gene_of)) < 2) {
  write_json_file(list(status = "skipped", reason = "Fewer than 2 genes with two or more expressed transcripts after filtering."),
                  file.path(out_dir, "dui_summary.json"))
  log_warn("Nothing to test for isoform usage.")
  quit(save = "no", status = 0)
}

# ---- DEXSeq ------------------------------------------------------------------------
mf <- model_frame(meta, params)
sd <- mf$col_data
covs_s <- make.names(params$covariates)
cond_s <- make.names(cond)
full <- as.formula(paste("~ sample + exon +", paste(c(paste0(covs_s, ":exon"), paste0(cond_s, ":exon")), collapse = " + ")))
reduced <- as.formula(paste("~ sample + exon", if (length(covs_s)) paste("+", paste0(covs_s, ":exon", collapse = " + ")) else ""))
bp <- if (params$threads > 1 && .Platform$OS.type == "unix") BiocParallel::MulticoreParam(params$threads) else BiocParallel::SerialParam()
dxd <- suppressMessages(DEXSeq::DEXSeqDataSet(countData = round(cnt), sampleData = sd, design = full,
                                              featureID = rownames(cnt), groupID = gene_of))
dxd <- DEXSeq::estimateSizeFactors(dxd)
dxd <- with_logged_warnings(DEXSeq::estimateDispersions(dxd, formula = full, quiet = TRUE, BPPARAM = bp))
dxd <- with_logged_warnings(DEXSeq::testForDEU(dxd, fullModel = full, reducedModel = reduced, BPPARAM = bp))
dxr <- DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)
gene_q <- DEXSeq::perGeneQValue(dxr)
res <- data.frame(transcript_id = dxr$featureID, gene_id = dxr$groupID, pvalue = dxr$pvalue, stringsAsFactors = FALSE)

# ---- Two-stage testing --------------------------------------------------------------
if (use_stager) {
  pconf <- matrix(ifelse(is.na(res$pvalue), 1, res$pvalue), ncol = 1, dimnames = list(res$transcript_id, "transcript"))
  st <- stageR::stageRTx(pScreen = gene_q, pConfirmation = pconf, pScreenAdjusted = TRUE,
                         tx2gene = data.frame(res$transcript_id, res$gene_id))
  st <- stageR::stageWiseAdjustment(st, method = "dtu", alpha = params$dui_alpha, allowNA = TRUE)
  adj <- suppressWarnings(stageR::getAdjustedPValues(st, order = FALSE, onlySignificantGenes = FALSE))
  res$gene_padj <- adj$gene[match(res$transcript_id, adj$txID)]
  res$transcript_padj <- adj$transcript[match(res$transcript_id, adj$txID)]
  method_txt <- sprintf("stageR %s two-stage (screen genes with DEXSeq per-gene q-values, confirm transcripts), overall FDR %g",
                        packageVersion("stageR"), params$dui_alpha)
} else {
  res$gene_padj <- unname(gene_q[res$gene_id])
  res$transcript_padj <- p.adjust(res$pvalue, "BH")
  method_txt <- "DEXSeq per-gene q-values and Benjamini-Hochberg per transcript (stageR not installed)"
}

# ---- Usage per group ----------------------------------------------------------------
prop <- share(cnt, gene_of)[res$transcript_id, , drop = FALSE]
g <- meta[[cond]]
res$usage_ref <- round(rowMeans(prop[, g == ref, drop = FALSE]), 4)
res$usage_test <- round(rowMeans(prop[, g == test, drop = FALSE]), 4)
res$delta_usage <- round(res$usage_test - res$usage_ref, 4)
names(res)[names(res) == "usage_ref"] <- paste0("usage_", ref)
names(res)[names(res) == "usage_test"] <- paste0("usage_", test)
res$gene_name <- tx$gene_name[match(res$transcript_id, tx$transcript_id)]
res$significant <- !is.na(res$gene_padj) & res$gene_padj < params$dui_alpha & !is.na(res$transcript_padj) &
  res$transcript_padj < params$dui_alpha & abs(res$delta_usage) >= params$min_delta_usage
res$direction <- ifelse(res$significant & res$delta_usage > 0, "up", ifelse(res$significant & res$delta_usage < 0, "down", "not significant"))
res <- res[order(res$gene_padj, res$transcript_padj, -abs(res$delta_usage)), c("transcript_id", "gene_id", "gene_name",
          paste0("usage_", ref), paste0("usage_", test), "delta_usage", "pvalue", "gene_padj", "transcript_padj", "significant", "direction")]
write_table(res, file.path(out_dir, "dui_transcripts_all.csv"))
write_table(res[res$significant, ], file.path(out_dir, "dui_transcripts_significant.csv"))

genes <- do.call(rbind, lapply(split(res, res$gene_id), function(d) {
  up <- d$transcript_id[d$direction == "up"]; dn <- d$transcript_id[d$direction == "down"]
  data.frame(gene_id = d$gene_id[1], gene_name = d$gene_name[1], n_transcripts_tested = nrow(d), gene_padj = d$gene_padj[1],
             n_transcripts_significant = sum(d$significant), max_abs_delta_usage = max(abs(d$delta_usage)),
             isoform_switch = length(up) > 0 && length(dn) > 0,
             transcripts_up = paste(up, collapse = ";"), transcripts_down = paste(dn, collapse = ";"),
             significant = any(d$significant), stringsAsFactors = FALSE)
}))
genes <- genes[order(!genes$significant, genes$gene_padj, -genes$max_abs_delta_usage), ]
write_table(genes, file.path(out_dir, "dui_genes.csv"))
rule <- sprintf("gene and transcript stage-wise adjusted p < %g and |delta usage| >= %g", params$dui_alpha, params$min_delta_usage)
log_info(sprintf("%d genes with significant isoform usage changes (%d isoform switches); %d transcripts; %s.",
                 sum(genes$significant), sum(genes$isoform_switch), sum(res$significant), rule))

# ---- Figures ----------------------------------------------------------------------------
dir_cols <- direction_colours(params)
v <- res[!is.na(res$transcript_padj), ]
v$neglog10 <- -log10(pmax(v$transcript_padj, 1e-300))
v$direction <- factor(v$direction, levels = names(dir_cols))
p <- ggplot(v[order(v$direction != "not significant"), ], aes(delta_usage, neglog10, colour = direction)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1) * params$min_delta_usage, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = dir_cols, drop = FALSE, name = NULL,
                      labels = c(up = sprintf("Larger share in %s", test), down = sprintf("Smaller share in %s", test), `not significant` = "Not significant")) +
  labs(title = sprintf("Differential isoform usage: %s vs %s", test, ref), subtitle = paste("Significant:", rule),
       x = sprintf("Change in transcript share of its gene (%s - %s)", test, ref), y = "-log10(adjusted p, transcript)") +
  theme_codebank() + theme(legend.position = "bottom")
save_plot(p, file.path(out_dir, "dui_volcano.png"), width = 8, height = 6)

top_genes <- head(genes$gene_id[genes$significant], params$n_usage_plots)
if (length(top_genes)) {
  rows <- res$gene_id %in% top_genes
  long <- do.call(rbind, lapply(which(rows), function(i) data.frame(
    gene = display_name(res$gene_name[i], res$gene_id[i]), transcript = res$transcript_id[i],
    sample = colnames(prop), group = g, usage = as.numeric(prop[res$transcript_id[i], ]), sig = res$significant[i])))
  long$group <- factor(long$group, levels = c(ref, test))
  long$gene <- factor(long$gene, levels = unique(display_name(genes$gene_name, genes$gene_id)[match(top_genes, genes$gene_id)]))
  means <- stats::aggregate(usage ~ gene + transcript + group + sig, data = long, FUN = mean)
  pal <- group_palette(c(ref, test), params)
  pu <- ggplot(means, aes(transcript, usage, fill = group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75, aes(alpha = sig)) +
    geom_point(data = long, aes(transcript, usage, group = group), position = position_dodge(width = 0.8), size = 0.9, show.legend = FALSE) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), labels = c(`TRUE` = "significant", `FALSE` = "not significant"), name = "Transcript") +
    scale_fill_manual(values = pal, name = NULL) + scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~gene, scales = "free_x", ncol = 3) +
    labs(title = sprintf("Transcript share of each gene, top %d genes with isoform usage changes", length(top_genes)),
         subtitle = "Bars: group means; points: samples", x = NULL, y = "Share of the gene's expression") +
    theme_codebank(10) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7), legend.position = "bottom")
  save_plot(pu, file.path(out_dir, "dui_top_genes_usage.png"), width = 11, height = 3 * ceiling(length(top_genes) / 3) + 1.5)
}

topg <- head(genes[genes$significant, ], 10)
write_json_file(list(
  method = sprintf("DEXSeq %s on dtuScaledTPM transcript counts (%s); %s",
                   as.character(packageVersion("DEXSeq")), paste(deparse(full), collapse = ""), method_txt),
  comparison = sprintf("%s vs %s", test, ref), samples_used = meta$sample_id,
  sign_convention = sprintf("delta_usage = share in %s - share in %s", test, ref),
  filter = sprintf("transcript >= %g reads and >= %g of its gene in >= %d samples; gene >= %g reads in every sample; >= 2 transcripts per gene",
                   params$min_feature_count, params$min_feature_prop, n_small, params$min_gene_count),
  significance_rule = rule, n_genes_tested = nrow(genes), n_transcripts_tested = nrow(res),
  n_genes_significant = sum(genes$significant), n_isoform_switches = sum(genes$isoform_switch),
  n_transcripts_significant = sum(res$significant),
  top_genes = if (nrow(topg)) lapply(seq_len(nrow(topg)), function(i) as.list(topg[i, c("gene_id", "gene_name", "gene_padj", "max_abs_delta_usage",
                                                                                      "isoform_switch", "transcripts_up", "transcripts_down")])) else list()
), file.path(out_dir, "dui_summary.json"))
write_step_provenance(run_dir, STEP, started,
                      params_used = params[c("min_feature_count", "min_feature_prop", "min_gene_count", "dui_alpha", "min_delta_usage", "threads")],
                      packages = c("DEXSeq", "stageR", "tximport", "ggplot2"), outputs = list.files(out_dir))
log_info("Step I03 done.")
