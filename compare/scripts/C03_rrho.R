#!/usr/bin/env Rscript
# =============================================================================
# C03_rrho.R  --  rank-rank hypergeometric overlap (RRHO) for every pair of studies
# =============================================================================
# Threshold-free: both studies' genes are ranked from most up to most down (by
# the metric chosen in C01) and the overlap of the two lists is tested at every
# pair of thresholds, in four quadrants (up/up, down/down, up/down, down/up; the
# stratified RRHO2 version). The map shows where the two rankings agree without
# picking a significance cutoff. P values are Benjamini-Yekutieli adjusted over
# the whole map (rrho_correction: BY).
# Outputs per pair: rrho_map.png, rrho_grid.csv, rrho_best_genes_<quadrant>.csv;
# rrho_summary.json with the strongest overlap of each quadrant.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "compare/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "compare", "R", "compare_utils.R"))
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

STEP <- "C03_rrho"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "Comparison step C03: rank-rank hypergeometric overlap."))
if (is_blank(opt[["run_dir"]])) fail("--run_dir is required.")
params <- load_compare_params(opt[["run_dir"]])
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
h <- read_harmonised(run_dir)
tables <- h$tables
correction <- params$rrho_correction %||% "BY"

QUAD_TEXT <- c(up_up = "up in both", down_down = "down in both", up_down = "up in first, down in second",
               down_up = "down in first, up in second")
pairs_out <- list()
for (pr in study_pairs(names(tables))) {
  a <- pr[1]; b <- pr[2]
  pdir <- step_dir(out_dir, pair_folder(a, b))
  ta <- tables[[a]][is.finite(tables[[a]]$rank_metric), ]
  tb <- tables[[b]][is.finite(tables[[b]]$rank_metric), ]
  va <- setNames(ta$rank_metric, ta$key); vb <- setNames(tb$rank_metric, tb$key)
  r <- tryCatch(rrho_map(va, vb, step = params$rrho_step, correction = correction), error = function(e) e)
  if (inherits(r, "error")) {
    log_warn(sprintf("RRHO %s vs %s skipped: %s", a, b, conditionMessage(r)))
    pairs_out[[length(pairs_out) + 1]] <- list(study_a = a, study_b = b, status = "skipped", reason = conditionMessage(r))
    next
  }
  g <- r$grid
  write_table(transform(g, pvalue = signif(pvalue, 4), padj = signif(padj, 4), neg_log10_padj = round(neg_log10_padj, 3)),
              file.path(pdir, "rrho_grid.csv"))
  gene_of <- setNames(ta$gene, ta$key)
  best <- lapply(names(r$best), function(q) {
    x <- r$best[[q]]
    if (is.null(x)) return(NULL)
    if (length(x$genes)) {
      bg <- data.frame(key = x$genes, gene = unname(gene_of[x$genes]),
                       rank_metric_a = round(unname(va[x$genes]), 3), rank_metric_b = round(unname(vb[x$genes]), 3))
      bg <- bg[order(-abs(bg$rank_metric_a) - abs(bg$rank_metric_b)), ]
      write_table(bg, file.path(pdir, paste0("rrho_best_genes_", q, ".csv")))
    }
    x$genes_file <- if (length(x$genes)) paste0("rrho_best_genes_", q, ".csv") else NULL
    x$top_genes <- head(unname(gene_of[x$genes][order(-abs(va[x$genes]) - abs(vb[x$genes]))]), 15)
    x$genes <- NULL
    x$meaning <- QUAD_TEXT[[q]]
    x
  })
  names(best) <- names(r$best)

  # Map: x = first study's ranking (most up on the left), y = second study's (most up at the bottom)
  cap <- max(g$neg_log10_padj, 1, na.rm = TRUE)
  p <- ggplot(g, aes(threshold_a, threshold_b, fill = neg_log10_padj)) + geom_raster() +
    scale_fill_gradientn(colours = c("#313695", "#4575B4", "#74ADD1", "#E0F3F8", "#FEE090", "#F46D43", "#A50026"),
                         limits = c(0, cap), name = if (correction == "BY") "-log10 adj. p\n(BY)" else "-log10 p") +
    geom_vline(xintercept = r$n_up_a, colour = "grey20", linewidth = 0.3, linetype = "dashed") +
    geom_hline(yintercept = r$n_up_b, colour = "grey20", linewidth = 0.3, linetype = "dashed") +
    annotate("text", x = r$n_up_a / 2, y = r$n_up_b * 0.75, label = "up in both", size = 3.2, colour = "white", fontface = "bold") +
    annotate("text", x = (r$n_up_a + r$n_genes) / 2, y = r$n_up_b + (r$n_genes - r$n_up_b) * 0.25, label = "down in both", size = 3.2, colour = "white", fontface = "bold") +
    scale_x_continuous(expand = c(0, 0), breaks = r$n_genes * c(0.04, 0.96), labels = c("most up", "most down")) +
    scale_y_continuous(expand = c(0, 0), breaks = r$n_genes * c(0.04, 0.96), labels = c("most up", "most down")) +
    labs(title = sprintf("Rank-rank overlap: %s vs %s", a, b),
         subtitle = sprintf("%d shared genes ranked by %s (step %d)\nDashed lines: up | down. Strongest: up in both %s (%d genes), down in both %s (%d genes)",
                            r$n_genes, h$info$rank_by, r$step,
                            format(best$up_up$padj, digits = 2), best$up_up$overlap, format(best$down_down$padj, digits = 2), best$down_down$overlap),
         x = paste(a, "ranking"), y = paste(b, "ranking")) +
    coord_fixed() + theme_codebank(10) + theme(panel.grid = element_blank())
  save_plot(p, file.path(pdir, "rrho_map.png"), width = 7.5, height = 6.8)
  concordant <- max(best$up_up$neg_log10_padj %||% 0, best$down_down$neg_log10_padj %||% 0)
  discordant <- max(best$up_down$neg_log10_padj %||% 0, best$down_up$neg_log10_padj %||% 0)
  pairs_out[[length(pairs_out) + 1]] <- list(
    study_a = a, study_b = b, status = "done", folder = basename(pdir), n_genes = r$n_genes, step = r$step,
    best = best, strongest_concordant = concordant, strongest_discordant = discordant,
    pattern = if (concordant < 1.3 && discordant < 1.3) "no significant overlap in any quadrant"
              else if (concordant >= discordant) "the studies agree (strongest overlap in the same-direction quadrants)"
              else "the studies disagree (strongest overlap in the opposite-direction quadrants)")
  log_info(sprintf("RRHO %s vs %s: up/up -log10 p %.1f (%d genes), down/down %.1f (%d genes).", a, b,
                   best$up_up$neg_log10_padj, best$up_up$overlap, best$down_down$neg_log10_padj, best$down_down$overlap))
}
write_json_file(list(method = sprintf("Stratified rank-rank hypergeometric overlap (RRHO2), %s correction", correction),
                     rank_by = h$info$rank_by, pairs = pairs_out), file.path(out_dir, "rrho_summary.json"))
write_step_provenance(run_dir, STEP, started, params_used = params[c("rrho_step", "rrho_correction")],
                      packages = c("ggplot2"), outputs = list.files(out_dir, recursive = TRUE))
