#!/usr/bin/env Rscript
# =============================================================================
# 05_de_plots.R  --  volcano, MA, heatmap and p-value histogram
# =============================================================================
# What it does
#   Draws the standard differential-expression figures from step 04's table.
#
#   volcano_plot.png       shrunken log2FC (x) vs -log10 p-value (y). Colours follow
#                          the significance rule (padj and |log2FC|). The dashed
#                          horizontal line is the raw p-value that corresponds to the
#                          padj cutoff. The top n_label_volcano genes (by padj) are labelled,
#                          or exactly the genes in volcano_label_genes when given.
#   ma_plot.png            mean normalised expression (x, log10) vs shrunken log2FC (y)
#   heatmap_top_degs.png   row-scaled (z-score) VST expression of the top
#                          n_heatmap_genes significant genes (or the top genes by
#                          padj if fewer are significant), samples clustered
#   pvalue_histogram.png   distribution of raw p-values (a check of the test:
#                          a flat histogram with a peak near 0 is the healthy shape)
#
# Outputs (in <run_dir>/05_de_plots/): the four PNGs and plots_summary.json
#
# Usage
#   Rscript scripts/05_de_plots.R --run_dir runs/my_run [--n_label_volcano 20 --n_heatmap_genes 40]
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "scripts/x"), ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(pheatmap)
  library(SummarizedExperiment)
})

STEP <- "05_de_plots"
started <- Sys.time()

option_list <- list(
  make_option("--run_dir", type = "character", help = "Run folder (required)"),
  make_option("--config", type = "character", help = "Parameter file (default: <run_dir>/params.yaml)"),
  make_option("--set", type = "character", help = "Extra overrides"),
  make_option("--n_label_volcano", type = "integer", help = "Genes labelled on the volcano plot"),
  make_option("--n_heatmap_genes", type = "integer", help = "Genes in the heatmap")
)
opt <- parse_args(OptionParser(option_list = option_list, description = "Step 05: DE figures."))
params <- load_params(opt)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
log_info("Step 05: DE plots")

de_path <- file.path(run_dir, "04_differential_expression", "de_results_all.csv")
if (!file.exists(de_path)) fail("Step 04 results not found: ", de_path)
de <- data.table::fread(de_path, data.table = FALSE, colClasses = list(character = c("gene_id", "gene_symbol", "entrez_id")))
de$label <- display_name(de$gene_symbol, de$gene_id)
ref <- params$reference_level; test <- params$test_level
pc <- params$padj_cutoff; fc <- params$log2fc_cutoff
dir_cols <- direction_colours(params)
n_up <- sum(de$direction == "up"); n_down <- sum(de$direction == "down")
title_cmp <- sprintf("%s vs %s", test, ref)
rule_txt <- sprintf("padj < %g and |log2FC| > %g", pc, fc)
has_repel <- requireNamespace("ggrepel", quietly = TRUE)

# ---- Volcano ----------------------------------------------------------------
v <- de[!is.na(de$pvalue), ]
min_p <- min(v$pvalue[v$pvalue > 0], na.rm = TRUE)
v$neglog10p <- -log10(pmax(v$pvalue, min_p))
v$direction <- factor(v$direction, levels = names(dir_cols))
v <- v[order(v$direction != "not significant"), ]   # draw significant points on top
p_line <- suppressWarnings(max(de$pvalue[!is.na(de$padj) & de$padj < pc], na.rm = TRUE))
# Labels: the genes the user named (volcano_label_genes), otherwise the top genes by padj
wanted <- unlist(params$volcano_label_genes %||% list())
if (length(wanted)) {
  hit <- toupper(v$label) %in% toupper(wanted) | toupper(v$gene_id) %in% toupper(wanted)
  to_label <- v[hit, ]
  not_found <- wanted[!toupper(wanted) %in% toupper(c(v$label, v$gene_id))]
  if (length(not_found)) log_warn("volcano_label_genes not among the tested genes: ", paste(not_found, collapse = ", "))
  label_note <- sprintf("Labelled: %d chosen genes.", nrow(to_label))
} else {
  not_found <- character(0)
  to_label <- head(de[de$significant, ], params$n_label_volcano)
  to_label <- v[v$gene_id %in% to_label$gene_id, ]
  label_note <- sprintf("Labelled: top %d by padj.", nrow(to_label))
}

p_vol <- ggplot(v, aes(log2FoldChange, neglog10p, colour = direction)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_vline(xintercept = c(-fc, fc), linetype = "dashed", colour = "grey40") +
  { if (is.finite(p_line)) geom_hline(yintercept = -log10(p_line), linetype = "dashed", colour = "grey40") } +
  scale_colour_manual(values = dir_cols, drop = FALSE,
                      labels = c(up = sprintf("Up in %s (%d)", test, n_up),
                                 down = sprintf("Down in %s (%d)", test, n_down),
                                 `not significant` = "Not significant"), name = NULL) +
  labs(title = sprintf("Volcano plot: %s", title_cmp),
       subtitle = sprintf("Significant: %s. %s", rule_txt, label_note),
       x = sprintf("log2 fold change (%s / %s, shrunken)", test, ref), y = "-log10(p-value)") +
  theme_codebank() + theme(legend.position = "bottom")
if (nrow(to_label)) {
  p_vol <- p_vol + if (has_repel) {
    ggrepel::geom_text_repel(data = to_label, aes(label = label), colour = "black", size = 3,
                             max.overlaps = 50, min.segment.length = 0, show.legend = FALSE)
  } else geom_text(data = to_label, aes(label = label), colour = "black", size = 3, vjust = -0.8, show.legend = FALSE)
}
save_plot(p_vol, file.path(out_dir, "volcano_plot.png"), width = 8, height = 6.5)

# ---- MA ---------------------------------------------------------------------
m <- de[de$baseMean > 0, ]
m$direction <- factor(m$direction, levels = names(dir_cols))
m <- m[order(m$direction != "not significant"), ]
p_ma <- ggplot(m, aes(log10(baseMean), log2FoldChange, colour = direction)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey30") +
  geom_hline(yintercept = c(-fc, fc), linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = dir_cols, drop = FALSE, name = NULL) +
  labs(title = sprintf("MA plot: %s", title_cmp), subtitle = sprintf("Significant: %s", rule_txt),
       x = "log10 mean normalised count", y = "log2 fold change (shrunken)") +
  theme_codebank() + theme(legend.position = "bottom")
save_plot(p_ma, file.path(out_dir, "ma_plot.png"), width = 7.5, height = 5.5)

# ---- Heatmap of top DEGs ----------------------------------------------------
vsd <- readRDS(file.path(run_dir, "03_exploratory", "vsd.rds"))
mat <- assay(vsd)
inp <- read_validated_inputs(run_dir)
meta <- inp$meta
grp <- condition_factor(meta, params)
sig_tab <- de[de$significant, ]
heat_note <- sprintf("top %d significant genes by padj", min(params$n_heatmap_genes, nrow(sig_tab)))
if (nrow(sig_tab) < 2) {
  sig_tab <- de[!is.na(de$padj), ]
  heat_note <- sprintf("fewer than 2 significant genes, so the top %d genes by padj are shown",
                       min(params$n_heatmap_genes, nrow(sig_tab)))
}
hm_genes <- head(sig_tab$gene_id, params$n_heatmap_genes)
hm_genes <- hm_genes[hm_genes %in% rownames(mat)]
heat_made <- FALSE
if (length(hm_genes) >= 2) {
  hm <- mat[hm_genes, , drop = FALSE]
  hm <- t(scale(t(hm)))
  hm[is.na(hm)] <- 0
  rownames(hm) <- de$label[match(hm_genes, de$gene_id)]
  ord <- order(grp)
  ann <- data.frame(row.names = colnames(hm), check.names = FALSE)
  ann[[params$condition]] <- grp
  for (cv in params$covariates) if (!cv %in% params$numeric_covariates) ann[[cv]] <- factor(meta[[cv]])
  lim <- max(2, ceiling(max(abs(hm))))
  pheatmap(hm[, ord, drop = FALSE], annotation_col = ann,
           annotation_colors = setNames(list(group_palette(levels(grp), params)), params$condition),
           color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
           breaks = seq(-lim, lim, length.out = 102),
           cluster_cols = TRUE, cluster_rows = TRUE, show_rownames = length(hm_genes) <= 80,
           fontsize_row = if (length(hm_genes) > 50) 6 else 8, border_color = NA,
           main = sprintf("%s\n%s; row z-scores of VST counts", title_cmp, heat_note),
           filename = file.path(out_dir, "heatmap_top_degs.png"),
           width = max(7, 0.3 * ncol(hm) + 4), height = max(6, 0.14 * length(hm_genes) + 3))
  heat_made <- TRUE
  log_info("Saved plot: heatmap_top_degs.png")
} else {
  log_warn("Not enough genes for a heatmap.")
}

# ---- P-value histogram ------------------------------------------------------
p_hist <- ggplot(de[!is.na(de$pvalue), ], aes(pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.025), fill = "#0072B2", colour = "white") +
  labs(title = "Distribution of p-values",
       subtitle = "Expected: flat, with a peak near 0 when there are real differences",
       x = "p-value", y = "Number of genes") +
  theme_codebank()
save_plot(p_hist, file.path(out_dir, "pvalue_histogram.png"), width = 6.5, height = 4.5)

write_json_file(list(plots = list.files(out_dir, pattern = "\\.png$"),
                     volcano_labelled_genes = to_label$label,
                     volcano_genes_not_found = not_found,
                     colours = as.list(dir_cols),
                     volcano_pvalue_line = if (is.finite(p_line)) signif(p_line, 4) else NULL,
                     heatmap = if (heat_made) list(n_genes = length(hm_genes), note = heat_note) else NULL),
                file.path(out_dir, "plots_summary.json"))
write_step_provenance(run_dir, STEP, started,
                      params_used = params[intersect(c("padj_cutoff", "log2fc_cutoff", "n_label_volcano", "n_heatmap_genes",
                                                       "volcano_label_genes", "up_colour", "down_colour",
                                                       "plot_palette", "group_colours"), names(params))],
                      packages = c("ggplot2", "ggrepel", "pheatmap"),
                      inputs = de_path, outputs = list.files(out_dir))
log_info("Step 05 done.")
