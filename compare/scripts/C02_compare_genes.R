#!/usr/bin/env Rscript
# =============================================================================
# C02_compare_genes.R  --  fold-change agreement and overlap of significant genes
# =============================================================================
# For every pair of studies, on the genes both studies tested:
#   - correlation of log2 fold changes (Pearson and Spearman; all shared genes,
#     and the genes significant in at least one study)
#   - overlap of significant genes by direction (up/up, down/down, up/down,
#     down/up), with one-sided Fisher's exact tests against the shared genes as
#     the universe, the Jaccard index and the share with the same direction
#   - a fold-change scatter plot and an overlap matrix figure
# Across all studies: a correlation heatmap, genes significant in two or more
# studies, and (3+ studies) the sizes of every intersection of significant genes.
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

STEP <- "C02_compare_genes"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "Comparison step C02: fold changes and significant genes."))
if (is_blank(opt[["run_dir"]])) fail("--run_dir is required.")
params <- load_compare_params(opt[["run_dir"]])
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
set.seed(params$seed)
h <- read_harmonised(run_dir)
tables <- h$tables
labels <- names(tables)
info <- h$info$studies
dcol <- direction_colours(params)
has_repel <- requireNamespace("ggrepel", quietly = TRUE)
warnings_out <- character(0)

cor_safe <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(list(r = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE))
  list(r = round(unname(ct$estimate), 4), p = signif(ct$p.value, 3), n = sum(ok))
}
fisher_greater <- function(in_a, in_b) {
  if (!any(in_a) || !any(in_b)) return(list(p = NA_real_, odds_ratio = NA_real_))
  ft <- stats::fisher.test(table(factor(in_a, c(TRUE, FALSE)), factor(in_b, c(TRUE, FALSE))), alternative = "greater")
  list(p = signif(ft$p.value, 3), odds_ratio = round(unname(ft$estimate), 3))
}
top_list <- function(d, n = 10) {
  if (!nrow(d)) return(list())
  d <- head(d[order(-pmin(abs(d$log2FC_a), abs(d$log2FC_b))), ], n)
  lapply(seq_len(nrow(d)), function(i) list(gene = d$gene[i], log2FC_a = round(d$log2FC_a[i], 2), log2FC_b = round(d$log2FC_b[i], 2)))
}

pair_results <- list()
for (pr in study_pairs(labels)) {
  a <- pr[1]; b <- pr[2]
  pdir <- step_dir(out_dir, pair_folder(a, b))
  ta <- tables[[a]]; tb <- tables[[b]]
  d <- merge(ta[, c("key", "gene", "log2FoldChange", "padj", "rank_metric", "significant", "direction")],
             tb[, c("key", "gene", "log2FoldChange", "padj", "rank_metric", "significant", "direction")],
             by = "key", suffixes = c("_a", "_b"))
  names(d) <- sub("^log2FoldChange_", "log2FC_", names(d))
  d$gene <- ifelse(is.na(d$gene_a) | d$gene_a == "", d$gene_b, d$gene_a)
  d$category <- ifelse(d$significant_a & d$significant_b,
                       ifelse(d$direction_a == d$direction_b, "both, same direction", "both, opposite directions"),
                       ifelse(d$significant_a, paste(a, "only"), ifelse(d$significant_b, paste(b, "only"), "neither")))
  n <- nrow(d)
  either <- d$significant_a | d$significant_b
  pear <- cor_safe(d$log2FC_a, d$log2FC_b, "pearson")
  spear <- cor_safe(d$log2FC_a, d$log2FC_b, "spearman")
  spear_rank <- cor_safe(d$rank_metric_a, d$rank_metric_b, "spearman")
  spear_sig <- cor_safe(d$log2FC_a[either], d$log2FC_b[either], "spearman")
  lv <- c("up", "down", "not significant")
  tab <- table(factor(d$direction_a, lv), factor(d$direction_b, lv))
  ov <- as.data.frame.matrix(tab)
  ov <- cbind(data.frame(direction_in_a = rownames(ov)), ov)
  names(ov) <- c(paste0("direction_in_", a), paste0(lv, "_in_", b))
  write_table(ov, file.path(pdir, "overlap_table.csv"))
  uu <- fisher_greater(d$direction_a == "up", d$direction_b == "up")
  dd <- fisher_greater(d$direction_a == "down", d$direction_b == "down")
  ud <- fisher_greater(d$direction_a == "up", d$direction_b == "down")
  du <- fisher_greater(d$direction_a == "down", d$direction_b == "up")
  n_both <- sum(d$significant_a & d$significant_b)
  n_same <- sum(d$category == "both, same direction")
  n_union <- sum(either)
  res <- list(
    study_a = a, study_b = b, folder = basename(pdir), n_shared_genes = n,
    n_significant_a = sum(d$significant_a), n_significant_b = sum(d$significant_b),
    n_significant_both = n_both, n_same_direction = n_same, n_opposite_direction = n_both - n_same,
    up_up = tab["up", "up"], down_down = tab["down", "down"], up_down = tab["up", "down"], down_up = tab["down", "up"],
    expected_up_up = round(sum(d$direction_a == "up") * sum(d$direction_b == "up") / max(n, 1), 1),
    expected_down_down = round(sum(d$direction_a == "down") * sum(d$direction_b == "down") / max(n, 1), 1),
    fisher_up_up = uu, fisher_down_down = dd, fisher_up_down = ud, fisher_down_up = du,
    jaccard = if (n_union) round(n_both / n_union, 4) else NA_real_,
    direction_agreement = if (n_both) round(n_same / n_both, 4) else NA_real_,
    pearson_log2FC = pear, spearman_log2FC = spear, spearman_rank_metric = spear_rank,
    spearman_log2FC_significant_either = spear_sig,
    top_shared_up = top_list(d[d$direction_a == "up" & d$direction_b == "up", ]),
    top_shared_down = top_list(d[d$direction_a == "down" & d$direction_b == "down", ]),
    top_opposite = top_list(d[d$category == "both, opposite directions", ])
  )
  if (is.finite(spear$r) && spear$r < -0.2)
    warnings_out <- c(warnings_out, sprintf("The fold changes of %s and %s are negatively correlated (Spearman %.2f). Check that both compare the groups the same way round (e.g. KO vs WT); if not, set flip: true for one of them.", a, b, spear$r))
  if (n_both > 0 && n_same / n_both < 0.5)
    warnings_out <- c(warnings_out, sprintf("Fewer than half of the genes significant in both %s and %s change in the same direction.", a, b))

  out <- d[order(d$category == "neither", -pmax(abs(d$rank_metric_a), abs(d$rank_metric_b))),
           c("key", "gene", "log2FC_a", "padj_a", "direction_a", "log2FC_b", "padj_b", "direction_b", "category")]
  names(out) <- c("key", "gene", paste0(c("log2FC_", "padj_", "direction_"), a)[1:3], paste0(c("log2FC_", "padj_", "direction_"), b), "category")
  write_table(out, file.path(pdir, "genes_compared.csv"))

  # ---- Figures ----
  cats <- c("both, same direction", "both, opposite directions", paste(a, "only"), paste(b, "only"), "neither")
  ccol <- setNames(c(dcol[["up"]], "#CC79A7", "#E69F00", dcol[["down"]], "grey78"), cats)
  d$category <- factor(d$category, cats)
  d <- d[order(d$category == "neither", decreasing = TRUE), ]
  same <- d[d$category == "both, same direction", ]
  lab <- head(same[order(-pmin(abs(same$log2FC_a), abs(same$log2FC_b))), ], params$n_label_genes)
  lim <- stats::quantile(abs(c(d$log2FC_a, d$log2FC_b)), 0.999, na.rm = TRUE)
  p <- ggplot(d, aes(log2FC_a, log2FC_b, colour = category)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.3) +
    geom_point(size = 0.9, alpha = 0.7) +
    scale_colour_manual(values = ccol, drop = TRUE, name = "Significant in") +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(title = sprintf("Fold changes: %s vs %s", a, b),
         subtitle = sprintf("%d shared genes; Spearman %.2f, Pearson %.2f; %d significant in both (%d same direction)",
                            n, spear$r, pear$r, n_both, n_same),
         x = sprintf("log2 fold change, %s (%s)", a, info$comparison[info$label == a]),
         y = sprintf("log2 fold change, %s (%s)", b, info$comparison[info$label == b])) +
    theme_codebank(11)
  if (nrow(lab) && has_repel) p <- p + ggrepel::geom_text_repel(data = lab, aes(label = gene), colour = "black", size = 3, max.overlaps = 30, show.legend = FALSE)
  save_plot(p, file.path(pdir, "fc_scatter.png"), width = 8, height = 6.5)

  m <- as.data.frame(as.table(tab)); names(m) <- c("a", "b", "n")
  m$label <- format(m$n, big.mark = ",")
  p <- ggplot(m, aes(b, a, fill = n)) + geom_tile(colour = "white") + geom_text(aes(label = label), size = 4) +
    scale_fill_gradient(low = "#F7FBFF", high = "#6BAED6", guide = "none") +
    scale_y_discrete(limits = rev(lv)) +
    labs(title = sprintf("Significant genes: %s vs %s", a, b),
         subtitle = sprintf("Up/up p = %s, down/down p = %s (one-sided Fisher's exact test)", format(uu$p, digits = 2), format(dd$p, digits = 2)),
         x = paste("In", b), y = paste("In", a)) + theme_codebank(11)
  save_plot(p, file.path(pdir, "overlap_matrix.png"), width = 6.5, height = 4.5)
  pair_results[[length(pair_results) + 1]] <- res
  log_info(sprintf("%s vs %s: %d shared genes, Spearman %.2f; %d significant in both (%d same direction).", a, b, n, spear$r, n_both, n_same))
}

# ---- Across all studies ---------------------------------------------------------------
k <- length(labels)
cm <- matrix(1, k, k, dimnames = list(labels, labels))
for (r in pair_results) { cm[r$study_a, r$study_b] <- cm[r$study_b, r$study_a] <- r$spearman_log2FC$r }
write_table(cbind(data.frame(study = labels), as.data.frame(cm, check.names = FALSE)), file.path(out_dir, "correlation_matrix.csv"))
cmd <- as.data.frame(as.table(cm)); names(cmd) <- c("a", "b", "r")
p <- ggplot(cmd, aes(a, b, fill = r)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", r)), size = 4) +
  scale_fill_gradient2(low = dcol[["down"]], mid = "white", high = dcol[["up"]], midpoint = 0, limits = c(-1, 1), name = "Spearman") +
  labs(title = "Correlation of log2 fold changes", subtitle = "Spearman, on the genes each pair of studies shares", x = NULL, y = NULL) +
  theme_codebank(11) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p, file.path(out_dir, "correlation_heatmap.png"), width = 2.2 + 1.1 * k, height = 1.6 + 0.9 * k)

all_keys <- Reduce(union, lapply(tables, `[[`, "key"))
dir_mat <- sapply(tables, function(t) t$direction[match(all_keys, t$key)])
dir_mat <- matrix(dir_mat, ncol = k, dimnames = list(all_keys, labels))
n_up <- rowSums(dir_mat == "up", na.rm = TRUE); n_down <- rowSums(dir_mat == "down", na.rm = TRUE)
n_tested <- rowSums(!is.na(dir_mat))
gene_of <- Reduce(function(x, t) ifelse(is.na(x), t$gene[match(all_keys, t$key)], x), tables, init = rep(NA_character_, length(all_keys)))
recurrent <- data.frame(key = all_keys, gene = gene_of, n_studies_tested = n_tested, n_up = n_up, n_down = n_down,
                        consistent = (n_up == 0 | n_down == 0), stringsAsFactors = FALSE)
for (l in labels) {
  t <- tables[[l]]; recurrent[[paste0("log2FC_", l)]] <- round(t$log2FoldChange[match(all_keys, t$key)], 3)
}
recurrent <- recurrent[(n_up + n_down) >= 2, ]
recurrent <- recurrent[order(-pmax(recurrent$n_up, recurrent$n_down), !recurrent$consistent, recurrent$gene), ]
write_table(recurrent, file.path(out_dir, "genes_significant_in_several_studies.csv"))

intersections <- NULL
if (k >= 3) {
  shared_all <- n_tested == k
  for (dirn in c("up", "down")) {
    sig <- dir_mat[shared_all, , drop = FALSE] == dirn
    combo <- apply(sig, 1, function(x) if (any(x)) paste(labels[x], collapse = " & ") else NA_character_)
    tt <- sort(table(combo[!is.na(combo)]), decreasing = TRUE)
    if (length(tt)) intersections <- rbind(intersections, data.frame(direction = dirn, studies = names(tt), n_studies = lengths(strsplit(names(tt), " & ")), n_genes = as.integer(tt)))
  }
  if (!is.null(intersections)) {
    write_table(intersections, file.path(out_dir, "intersections.csv"))
    top <- do.call(rbind, lapply(split(intersections, intersections$direction), function(x) head(x, 12)))
    top$studies <- factor(top$studies, levels = rev(unique(top$studies[order(top$n_studies, top$n_genes)])))
    p <- ggplot(top, aes(n_genes, studies, fill = direction)) + geom_col() + facet_wrap(~direction, scales = "free_y") +
      scale_fill_manual(values = c(up = dcol[["up"]], down = dcol[["down"]]), guide = "none") +
      labs(title = "Significant genes shared by the studies", subtitle = sprintf("Genes tested in all %d studies; each gene is counted once, in the exact set of studies where it is significant", k),
           x = "Genes", y = NULL) + theme_codebank(10)
    save_plot(p, file.path(out_dir, "intersections.png"), width = 10, height = max(4, 0.35 * max(table(top$direction)) + 2))
  }
}

write_json_file(list(
  significance_rule = h$info$significance_rule, match_method = h$info$match_method,
  n_genes_in_all_studies = h$info$n_genes_in_all_studies,
  pairs = pair_results,
  n_genes_significant_in_several = nrow(recurrent), n_consistent_in_several = sum(recurrent$consistent),
  intersections = intersections, warnings = as.list(warnings_out)
), file.path(out_dir, "gene_comparison.json"))
for (w in warnings_out) log_warn(w)
write_step_provenance(run_dir, STEP, started, params_used = params[c("n_label_genes", "up_colour", "down_colour")],
                      packages = c("ggplot2", "ggrepel"), outputs = list.files(out_dir, recursive = TRUE))
