#!/usr/bin/env Rscript
# =============================================================================
# C04_compare_pathways.R  --  compare the studies' pathway results (GSEA)
# =============================================================================
# Uses each code-bank run's own GSEA results (step 06, fgsea on the Wald
# statistic). For a study given as a results table, GSEA is run here the same
# way (fgsea, same gene sets and size limits, ranked by the metric chosen in
# C01). Pathways are matched by name within each collection (MSigDB names are
# the same for human and mouse).
# Per collection: NES of every pathway in every study (table and heatmap of the
# top pathways); per pair: NES correlation, pathways significant in both with the
# same / opposite direction, NES scatter plot.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "compare/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "R", "genesets.R"))
source(file.path(REPO_ROOT, "compare", "R", "compare_utils.R"))
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

STEP <- "C04_compare_pathways"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "Comparison step C04: pathways across studies."))
if (is_blank(opt[["run_dir"]])) fail("--run_dir is required.")
params <- load_compare_params(opt[["run_dir"]])
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
set.seed(params$seed)
h <- read_harmonised(run_dir)
tables <- h$tables
info <- h$info$studies
labels <- names(tables)
cut <- params$pathway_padj_cutoff
dcol <- direction_colours(params)

finish_skipped <- function(reason) {
  log_warn("Pathway comparison skipped: ", reason)
  write_json_file(list(status = "skipped", reason = reason), file.path(out_dir, "pathway_comparison.json"))
  write_step_provenance(run_dir, STEP, started, params_used = params[c("compare_pathways", "gene_set_collections")],
                        outputs = list.files(out_dir), notes = list(status = "skipped", reason = reason))
  quit(save = "no", status = 0)
}
if (!isTRUE(params$compare_pathways)) finish_skipped("compare_pathways is false.")

folder_of <- function(key) make.names(gsub(":", "_", key))
wanted <- params$gene_set_collections
gsea <- list()           # gsea[[study]][[collection key]] = data.frame(pathway, NES, padj, set_size)
notes <- list()
for (l in labels) {
  s <- info[info$label == l, ]
  gsea[[l]] <- list()
  if (s$source == "run") {
    # every collection the run tested (MSigDB keys are stored as folder names: GO:BP -> GO_BP)
    for (f in list.files(file.path(s$input_dir, "06_enrichment"), pattern = "^gsea\\.csv$", recursive = TRUE, full.names = TRUE)) {
      g <- data.table::fread(f, data.table = FALSE)
      if (nrow(g)) gsea[[l]][[basename(dirname(f))]] <- g[, c("pathway", "NES", "padj", "set_size")]
    }
  } else {
    if (!requireNamespace("fgsea", quietly = TRUE)) { notes[[l]] <- "fgsea is not installed"; next }
    src <- params$gene_set_source %||% "msigdbr"
    use_msig <- src %in% c("msigdbr", "both") && s$organism %in% c("human", "mouse")
    use_gmt <- src %in% c("gmt", "both") && length(params$gmt_files) > 0
    if (!use_msig && !use_gmt) { notes[[l]] <- sprintf("no gene sets for organism '%s' (set organism, or give gmt_files)", s$organism); next }
    gs <- load_gene_set_collections(use_msig, use_gmt, wanted, params$gmt_files, s$organism)
    t <- tables[[l]]
    t <- t[is.finite(t$rank_metric) & !is.na(t$gene_symbol) & t$gene_symbol != "", ]
    t <- t[order(-abs(t$rank_metric)), ]; t <- t[!duplicated(t$gene_symbol), ]
    ranks <- setNames(t$rank_metric, t$gene_symbol)
    for (ckey in names(gs$collections)) {
      key <- folder_of(ckey)
      r <- suppressWarnings(fgsea::fgsea(pathways = gs$collections[[ckey]]$sets, stats = ranks,
                                         minSize = params$gene_set_min_size, maxSize = params$gene_set_max_size, eps = 0))
      r <- as.data.frame(r)[, c("pathway", "NES", "padj", "size")]
      names(r)[4] <- "set_size"
      gsea[[l]][[key]] <- r
      dir.create(file.path(out_dir, "gsea_run_here"), showWarnings = FALSE)
      write_table(r[order(r$padj), ], file.path(out_dir, "gsea_run_here", sprintf("%s_%s.csv", safe_label(l), folder_of(key))))
    }
    notes[[l]] <- sprintf("GSEA run here (fgsea %s, ranked by %s)", as.character(utils::packageVersion("fgsea")), h$info$rank_by)
    log_info(sprintf("%s: GSEA run on %d ranked genes for %s.", l, length(ranks), paste(names(gs$collections), collapse = ", ")))
  }
}
all_labels <- labels
labels <- all_labels[vapply(gsea, length, 0L) > 0]        # studies without any GSEA result are left out here
left_out <- setdiff(all_labels, labels)
if (length(labels) < 2) finish_skipped(sprintf("fewer than two studies have GSEA results (%s).",
                                               paste(sprintf("%s: %s", all_labels, vapply(all_labels, function(l) notes[[l]] %||% "none", "")), collapse = "; ")))
if (length(left_out)) log_warn("Left out of the pathway comparison (no GSEA results): ", paste(left_out, collapse = ", "))
common <- Reduce(intersect, lapply(gsea[labels], names))
order_pref <- vapply(wanted, folder_of, "")
common <- c(intersect(order_pref, common), setdiff(common, order_pref))
if (!length(common)) finish_skipped(sprintf("no gene-set collection has GSEA results in every study (%s).",
                                            paste(sprintf("%s: %s", labels, vapply(gsea[labels], function(x) if (length(x)) paste(names(x), collapse = ", ") else "none", "")), collapse = "; ")))

coll_out <- list()
pair_out <- list()
for (key in common) {
  cdir <- step_dir(out_dir, folder_of(key))
  all_pw <- Reduce(union, lapply(labels, function(l) gsea[[l]][[key]]$pathway))
  wide <- data.frame(pathway = all_pw, stringsAsFactors = FALSE)
  for (l in labels) {
    g <- gsea[[l]][[key]]; m <- match(all_pw, g$pathway)
    wide[[paste0("NES_", l)]] <- round(g$NES[m], 3)
    wide[[paste0("padj_", l)]] <- signif(g$padj[m], 3)
  }
  nes <- as.matrix(wide[, paste0("NES_", labels), drop = FALSE]); pad <- as.matrix(wide[, paste0("padj_", labels), drop = FALSE])
  sig <- !is.na(pad) & pad < cut                              # same row order as wide (before sorting)
  wide$n_significant_up <- rowSums(sig & nes > 0); wide$n_significant_down <- rowSums(sig & nes < 0)
  wide$consistent <- wide$n_significant_up == 0 | wide$n_significant_down == 0
  wide$min_padj <- apply(pad, 1, function(x) suppressWarnings(min(x, na.rm = TRUE)))
  wide <- wide[order(-(wide$n_significant_up + wide$n_significant_down), wide$min_padj), ]
  write_table(wide, file.path(cdir, "pathways_all_studies.csv"))

  top <- head(wide[(wide$n_significant_up + wide$n_significant_down) > 0, ], params$n_heatmap_pathways)
  if (nrow(top)) {
    hm <- do.call(rbind, lapply(labels, function(l) data.frame(pathway = top$pathway, study = l, NES = top[[paste0("NES_", l)]],
                                                               star = ifelse(!is.na(top[[paste0("padj_", l)]]) & top[[paste0("padj_", l)]] < cut, "*", ""))))
    hm$pathway <- factor(pretty_path(hm$pathway), levels = rev(unique(pretty_path(top$pathway))))
    hm$study <- factor(hm$study, labels)
    lim <- max(abs(hm$NES), 1, na.rm = TRUE)
    p <- ggplot(hm, aes(study, pathway, fill = NES)) + geom_tile(colour = "white") + geom_text(aes(label = star), size = 4) +
      scale_fill_gradient2(low = dcol[["down"]], mid = "white", high = dcol[["up"]], midpoint = 0, limits = c(-lim, lim), na.value = "grey90") +
      labs(title = sprintf("%s: pathways across studies", key),
           subtitle = sprintf("NES of the top %d pathways significant in at least one study; * padj < %g", nrow(top), cut), x = NULL, y = NULL) +
      theme_codebank(9) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
    save_plot(p, file.path(cdir, "nes_heatmap.png"), width = 7 + 0.6 * length(labels), height = max(4, 0.28 * nrow(top) + 2))
  }
  n_sig_studies <- wide$n_significant_up + wide$n_significant_down
  coll_out[[key]] <- list(folder = basename(cdir), n_pathways = nrow(wide),
                          n_significant_in_every_study = sum(n_sig_studies == length(labels)),
                          n_significant_in_two_or_more = sum(n_sig_studies >= 2),
                          n_consistent_in_two_or_more = sum(n_sig_studies >= 2 & wide$consistent),
                          heatmap = if (nrow(top)) "nes_heatmap.png" else NULL)

  for (pr in study_pairs(labels)) {
    a <- pr[1]; b <- pr[2]
    d <- data.frame(pathway = wide$pathway, NES_a = wide[[paste0("NES_", a)]], padj_a = wide[[paste0("padj_", a)]],
                    NES_b = wide[[paste0("NES_", b)]], padj_b = wide[[paste0("padj_", b)]])
    d <- d[!is.na(d$NES_a) & !is.na(d$NES_b), ]
    sa <- !is.na(d$padj_a) & d$padj_a < cut; sb <- !is.na(d$padj_b) & d$padj_b < cut
    d$category <- ifelse(sa & sb, ifelse(sign(d$NES_a) == sign(d$NES_b), "both, same direction", "both, opposite directions"),
                         ifelse(sa, paste(a, "only"), ifelse(sb, paste(b, "only"), "neither")))
    rr <- if (nrow(d) >= 10) suppressWarnings(stats::cor.test(d$NES_a, d$NES_b, method = "spearman", exact = FALSE)) else NULL
    pdir <- step_dir(cdir, pair_folder(a, b))
    out <- d[order(d$category == "neither", -pmin(abs(d$NES_a), abs(d$NES_b))), ]
    names(out) <- c("pathway", paste0("NES_", a), paste0("padj_", a), paste0("NES_", b), paste0("padj_", b), "category")
    write_table(out, file.path(pdir, "pathways_compared.csv"))
    cats <- c("both, same direction", "both, opposite directions", paste(a, "only"), paste(b, "only"), "neither")
    ccol <- setNames(c(dcol[["up"]], "#CC79A7", "#E69F00", dcol[["down"]], "grey78"), cats)
    d$category <- factor(d$category, cats)
    same <- d[d$category == "both, same direction", ]
    lab <- head(same[order(-pmin(abs(same$NES_a), abs(same$NES_b))), ], 8)
    p <- ggplot(d[order(d$category == "neither", decreasing = TRUE), ], aes(NES_a, NES_b, colour = category)) +
      geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
      geom_point(size = 1.4, alpha = 0.8) + scale_colour_manual(values = ccol, drop = TRUE, name = sprintf("padj < %g in", cut)) +
      labs(title = sprintf("%s pathways: %s vs %s", key, a, b),
           subtitle = sprintf("%d pathways in both; Spearman %s; %d significant in both (%d same direction)", nrow(d),
                              if (is.null(rr)) "NA" else sprintf("%.2f", rr$estimate), sum(sa & sb), sum(d$category == "both, same direction")),
           x = paste("NES,", a), y = paste("NES,", b)) + theme_codebank(10)
    if (nrow(lab) && requireNamespace("ggrepel", quietly = TRUE))
      p <- p + ggrepel::geom_text_repel(data = lab, aes(label = pretty_path(pathway, 35)), colour = "black", size = 2.6, max.overlaps = 20, show.legend = FALSE)
    save_plot(p, file.path(pdir, "nes_scatter.png"), width = 8.5, height = 6.5)
    top_same <- head(same[order(-pmin(abs(same$NES_a), abs(same$NES_b))), ], 8)
    pair_out[[length(pair_out) + 1]] <- list(
      collection = key, study_a = a, study_b = b, folder = file.path(basename(cdir), basename(pdir)),
      n_pathways = nrow(d), spearman_NES = if (is.null(rr)) NA_real_ else round(unname(rr$estimate), 4),
      spearman_p = if (is.null(rr)) NA_real_ else signif(rr$p.value, 3),
      n_significant_a = sum(sa), n_significant_b = sum(sb), n_significant_both = sum(sa & sb),
      n_same_direction = sum(d$category == "both, same direction"), n_opposite_direction = sum(d$category == "both, opposite directions"),
      top_shared = lapply(seq_len(nrow(top_same)), function(i) list(pathway = top_same$pathway[i], NES_a = top_same$NES_a[i], NES_b = top_same$NES_b[i])))
  }
  log_info(sprintf("%s: %d pathways; %d significant in two or more studies.", key, nrow(wide), coll_out[[key]]$n_significant_in_two_or_more))
}
write_json_file(list(status = "done", padj_cutoff = cut, collections_compared = common, studies_compared = labels,
                     studies_left_out = as.list(left_out),
                     collections_not_in_every_study = setdiff(unique(unlist(lapply(gsea, names))), common),
                     gsea_source = lapply(setNames(all_labels, all_labels), function(l) notes[[l]] %||% if (length(gsea[[l]])) "the run's own GSEA (step 06)" else "no GSEA results"),
                     collections = coll_out, pairs = pair_out), file.path(out_dir, "pathway_comparison.json"))
write_step_provenance(run_dir, STEP, started, params_used = params[c("compare_pathways", "pathway_padj_cutoff", "gene_set_collections",
                                                                   "gene_set_min_size", "gene_set_max_size")],
                      packages = c("fgsea", "msigdbr", "ggplot2"), outputs = list.files(out_dir, recursive = TRUE))
