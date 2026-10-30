#!/usr/bin/env Rscript
# =============================================================================
# S02_splicing_results.R  --  significant splicing events, tables and figures
# =============================================================================
# Statistics are rMATS's (likelihood-ratio test on inclusion levels, FDR by
# Benjamini-Hochberg within each event type). This step applies the run's rule:
#   significant = FDR <= fdr_cutoff and |dPSI| >= delta_psi_cutoff and
#                 mean reads per sample >= min_reads_per_sample in BOTH groups
# dPSI = PSI(group 1) - PSI(group 2): positive = more inclusion in group 1 ("up").
#
# Outputs (in <run_dir>/S02_splicing_results/):
#   splicing_events_all.csv          every event with coverage flag, significance, direction
#   splicing_events_significant.csv  significant events, strongest first
#   splicing_summary.json            counts per event type and direction, top events
#   event_counts.png, dpsi_volcano.png, psi_scatter.png, dpsi_distribution.png, top_events_heatmap.png
# With de_run (a main-pipeline run folder for the same comparison):
#   splicing_vs_expression.csv / .png   which genes change splicing, expression, or both
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "splicing/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "splicing", "R", "splicing_utils.R"))
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

STEP <- "S02_splicing_results"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "Splicing step S02: significant events and figures."))
if (is_blank(opt$run_dir)) fail("--run_dir is required.")
params <- load_splicing_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
log_info("Step S02: significant splicing events")

ev <- data.table::fread(file.path(run_dir, "S01_validate_rmats", "events_validated.csv"), data.table = FALSE,
                        colClasses = list(character = c("gene_id", "gene_symbol", "chr", "psi_values_1", "psi_values_2", "coordinates")))
g1 <- params$group_1_name; g2 <- params$group_2_name
fdr_cut <- params$fdr_cutoff; dpsi_cut <- params$delta_psi_cutoff; min_reads <- params$min_reads_per_sample

ev$passes_coverage <- !is.na(ev$mean_reads_1) & !is.na(ev$mean_reads_2) &
  ev$mean_reads_1 >= min_reads & ev$mean_reads_2 >= min_reads
ev$significant <- ev$passes_coverage & !is.na(ev$fdr) & !is.na(ev$delta_psi) &
  ev$fdr <= fdr_cut & abs(ev$delta_psi) >= dpsi_cut
ev$direction <- ifelse(ev$significant & ev$delta_psi > 0, "up", ifelse(ev$significant & ev$delta_psi < 0, "down", "not significant"))
ev$note <- ifelse(!ev$passes_coverage, sprintf("fewer than %g reads per sample in a group", min_reads),
                  ifelse(is.na(ev$fdr), "no test result", ""))
ev$label <- event_label(ev)
ev <- ev[order(ev$fdr, -abs(ev$delta_psi)), ]
rule <- sprintf("FDR <= %g and |dPSI| >= %g, with >= %g reads per sample in both groups", fdr_cut, dpsi_cut, min_reads)

write_table(ev, file.path(out_dir, "splicing_events_all.csv"))
sig <- ev[ev$significant, ]
write_table(sig, file.path(out_dir, "splicing_events_significant.csv"))

types <- params$event_types
count_rows <- lapply(types, function(t) {
  e <- ev[ev$event_type == t, ]
  data.frame(event_type = t, label = EVENT_LABELS[[t]], tested = nrow(e), passes_coverage = sum(e$passes_coverage),
             significant = sum(e$significant), up = sum(e$direction == "up"), down = sum(e$direction == "down"))
})
counts <- do.call(rbind, count_rows)
log_info(sprintf("%d of %d events significant (%s).", nrow(sig), nrow(ev), rule))

# ---- Figures ---------------------------------------------------------------------
dir_cols <- direction_colours(params)
dir_labels <- c(up = sprintf("More inclusion in %s", g1), down = sprintf("Less inclusion in %s", g1),
                `not significant` = "Not significant")

cl <- rbind(data.frame(event_type = counts$event_type, direction = "up", n = counts$up),
            data.frame(event_type = counts$event_type, direction = "down", n = counts$down))
cl$event_type <- factor(cl$event_type, levels = types)
cl$direction <- factor(cl$direction, levels = c("up", "down"))
p_counts <- ggplot(cl, aes(event_type, n, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = dir_cols[c("up", "down")], labels = dir_labels[c("up", "down")], name = NULL) +
  labs(title = sprintf("Significant splicing events: %s vs %s", g1, g2), subtitle = rule,
       x = NULL, y = "Events") +
  theme_codebank() + theme(legend.position = "bottom")
save_plot(p_counts, file.path(out_dir, "event_counts.png"), width = 8, height = 5)

vol <- ev[ev$passes_coverage & !is.na(ev$fdr) & !is.na(ev$delta_psi), ]
if (nrow(vol)) {
  min_fdr <- min(vol$fdr[vol$fdr > 0], 1e-300, na.rm = TRUE)
  vol$neglog10fdr <- -log10(pmax(vol$fdr, min_fdr))
  vol$direction <- factor(vol$direction, levels = names(dir_cols))
  vol$event_type <- factor(vol$event_type, levels = types)
  vol <- vol[order(vol$direction != "not significant"), ]
  lab <- head(vol[vol$significant, ][order(-vol$neglog10fdr[vol$significant], -abs(vol$delta_psi[vol$significant])), ],
              params$n_label_events)
  p_vol <- ggplot(vol, aes(delta_psi, neglog10fdr, colour = direction)) +
    geom_point(size = 1.1, alpha = 0.7) +
    geom_vline(xintercept = c(-dpsi_cut, dpsi_cut), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = dir_cols, labels = dir_labels, drop = FALSE, name = NULL) +
    facet_wrap(~event_type, nrow = 1) +
    labs(title = sprintf("Change in exon inclusion: %s vs %s", g1, g2),
         subtitle = sprintf("Events with enough reads; significant: %s", rule),
         x = sprintf("dPSI (%s - %s)", g1, g2), y = "-log10(FDR)") +
    theme_codebank() + theme(legend.position = "bottom")
  if (nrow(lab)) {
    lab$event_type <- factor(lab$event_type, levels = types)
    lab$gene <- display_name(lab$gene_symbol, lab$gene_id)
    p_vol <- p_vol + if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(data = lab, aes(label = gene), colour = "black", size = 2.8, max.overlaps = 30,
                               min.segment.length = 0, show.legend = FALSE)
    } else geom_text(data = lab, aes(label = gene), colour = "black", size = 2.8, vjust = -0.8, show.legend = FALSE)
  }
  save_plot(p_vol, file.path(out_dir, "dpsi_volcano.png"), width = max(8, 2.6 * length(types)), height = 5.5)

  p_sc <- ggplot(vol, aes(psi_2, psi_1, colour = direction)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed") +
    geom_point(size = 1, alpha = 0.7) +
    scale_colour_manual(values = dir_cols, labels = dir_labels, drop = FALSE, name = NULL) +
    facet_wrap(~event_type, nrow = 1) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = "Percent spliced in (PSI) per event, group means", x = sprintf("PSI in %s", g2), y = sprintf("PSI in %s", g1)) +
    theme_codebank() + theme(legend.position = "bottom")
  save_plot(p_sc, file.path(out_dir, "psi_scatter.png"), width = max(8, 2.6 * length(types)), height = 4.2)
}

# dPSI distribution per event type (all events with enough reads; significant ones on top)
if (nrow(vol)) {
  p_dist <- ggplot(vol, aes(event_type, delta_psi)) +
    geom_violin(fill = "grey90", colour = "grey60", scale = "width") +
    geom_jitter(data = vol[vol$significant, ], aes(colour = direction), width = 0.15, size = 1.2, alpha = 0.8) +
    geom_hline(yintercept = c(-dpsi_cut, 0, dpsi_cut), linetype = c("dashed", "solid", "dashed"), colour = "grey40") +
    scale_colour_manual(values = dir_cols, labels = dir_labels, drop = FALSE, name = NULL) +
    labs(title = "Change in inclusion by event type", subtitle = "Violins: all events with enough reads; points: significant events",
         x = NULL, y = sprintf("dPSI (%s - %s)", g1, g2)) +
    theme_codebank() + theme(legend.position = "bottom")
  save_plot(p_dist, file.path(out_dir, "dpsi_distribution.png"), width = 7.5, height = 5.5)
}

# ---- Splicing vs expression (optional: a main-pipeline run of the same comparison) -----
integration <- NULL
if (!is_blank(params$de_run)) {
  de_path <- file.path(resolve_path(params$de_run) %||% params$de_run, "04_differential_expression", "de_results_all.csv")
  if (!file.exists(de_path)) {
    log_warn("de_run has no 04_differential_expression/de_results_all.csv; splicing vs expression skipped.")
  } else {
    de <- data.table::fread(de_path, data.table = FALSE, colClasses = list(character = c("gene_id", "gene_symbol")))
    strip <- function(x) sub("\\.\\d+$", "", x)
    tested <- ev[ev$passes_coverage, ]
    by_gene <- do.call(rbind, lapply(split(tested, strip(tested$gene_id)), function(d) {
      s <- d[d$significant, ]
      top <- if (nrow(s)) s[which.max(abs(s$delta_psi)), ] else d[which.max(abs(d$delta_psi)), ]
      data.frame(gene_id = strip(d$gene_id[1]), gene_symbol = d$gene_symbol[1], n_events_tested = nrow(d),
                 n_events_significant = nrow(s), strongest_event = top$event_id, strongest_delta_psi = round(top$delta_psi, 3),
                 splicing_significant = nrow(s) > 0, stringsAsFactors = FALSE)
    }))
    m <- match(by_gene$gene_id, strip(de$gene_id))
    by_gene$de_log2FoldChange <- round(de$log2FoldChange[m], 3)
    by_gene$de_padj <- signif(de$padj[m], 3)
    by_gene$expression_significant <- !is.na(m) & de$significant[m] %in% c(TRUE, "TRUE", "True")
    by_gene$category <- ifelse(by_gene$splicing_significant & by_gene$expression_significant, "splicing and expression",
                        ifelse(by_gene$splicing_significant, "splicing only",
                        ifelse(by_gene$expression_significant, "expression only", "neither")))
    by_gene$category[is.na(m)] <- ifelse(by_gene$splicing_significant[is.na(m)], "splicing only (gene not in DE results)", "not in DE results")
    by_gene <- by_gene[order(by_gene$category == "neither", -by_gene$splicing_significant, -abs(by_gene$strongest_delta_psi)), ]
    write_table(by_gene, file.path(out_dir, "splicing_vs_expression.csv"))
    cats <- c("splicing and expression" = "#CC79A7", "splicing only" = "#009E73", "expression only" = "#E69F00", "neither" = "grey80")
    pd <- by_gene[by_gene$category %in% names(cats) & !is.na(by_gene$de_log2FoldChange), ]
    pd$category <- factor(pd$category, levels = names(cats))
    p_int <- ggplot(pd[order(pd$category == "neither", decreasing = TRUE), ], aes(de_log2FoldChange, strongest_delta_psi, colour = category)) +
      geom_hline(yintercept = 0, colour = "grey50") + geom_vline(xintercept = 0, colour = "grey50") +
      geom_point(size = 1.3, alpha = 0.8) + scale_colour_manual(values = cats, drop = FALSE, name = NULL) +
      labs(title = "Splicing changes vs expression changes per gene",
           subtitle = "x: gene log2 fold change (main pipeline); y: dPSI of the gene's strongest splicing event",
           x = sprintf("log2 fold change (%s / %s)", g1, g2), y = "dPSI of the strongest event") +
      theme_codebank() + theme(legend.position = "bottom")
    save_plot(p_int, file.path(out_dir, "splicing_vs_expression.png"), width = 7.5, height = 6)
    integration <- list(de_run = basename(dirname(dirname(de_path))), genes = as.list(table(by_gene$category)),
                        note = "Genes with splicing changes but no expression change point to regulation of isoforms rather than of the gene's level.")
    log_info("Splicing vs expression: ", paste(names(table(by_gene$category)), table(by_gene$category), sep = " ", collapse = "; "))
  }
}

heat_made <- FALSE
if (nrow(sig) >= 2 && requireNamespace("pheatmap", quietly = TRUE)) {
  top <- head(sig[order(sig$fdr, -abs(sig$delta_psi)), ], params$n_heatmap_events)
  m1 <- do.call(rbind, split_numbers(top$psi_values_1))
  m2 <- do.call(rbind, split_numbers(top$psi_values_2))
  mat <- cbind(m1, m2)
  if (ncol(mat) == length(params$group_1_samples) + length(params$group_2_samples)) {
    colnames(mat) <- c(params$group_1_samples, params$group_2_samples)
    rownames(mat) <- make.unique(sprintf("%s %s", top$event_type, display_name(top$gene_symbol, top$gene_id)))
    ann <- data.frame(group = rep(c(g1, g2), c(ncol(m1), ncol(m2))), row.names = colnames(mat))
    pal <- group_palette(c(g2, g1), params)
    pheatmap::pheatmap(mat, cluster_rows = nrow(mat) > 2, cluster_cols = FALSE, na_col = "grey85",
                       color = colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(100), breaks = seq(0, 1, length.out = 101),
                       annotation_col = ann, annotation_colors = list(group = pal), border_color = NA,
                       main = sprintf("PSI per sample, top %d significant events", nrow(mat)),
                       filename = file.path(out_dir, "top_events_heatmap.png"),
                       width = max(6, 0.4 * ncol(mat) + 4), height = max(5, 0.16 * nrow(mat) + 2.5))
    heat_made <- TRUE
    log_info("Saved plot: top_events_heatmap.png")
  }
}

# ---- Summary ------------------------------------------------------------------------
top_list <- function(d, n = 5) {
  d <- head(d, n)
  if (!nrow(d)) return(list())
  lapply(seq_len(nrow(d)), function(i) list(event_id = d$event_id[i], gene = display_name(d$gene_symbol[i], d$gene_id[i]),
                                           event_type = d$event_type[i], location = sprintf("%s:%d-%d (%s)", d$chr[i], d$event_start[i], d$event_end[i], d$strand[i]),
                                           delta_psi = round(d$delta_psi[i], 3), psi_1 = d$psi_1[i], psi_2 = d$psi_2[i],
                                           fdr = signif(d$fdr[i], 3)))
}
summary <- list(
  comparison = sprintf("%s vs %s", g1, g2),
  sign_convention = sprintf("dPSI = PSI(%s) - PSI(%s); positive = more inclusion of the alternative exon/region in %s", g1, g2, g1),
  significance_rule = rule,
  count_type = params$count_type,
  n_events_tested = nrow(ev), n_events_passing_coverage = sum(ev$passes_coverage),
  n_significant = nrow(sig), n_up = sum(sig$direction == "up"), n_down = sum(sig$direction == "down"),
  n_genes_significant = length(unique(na.omit(display_name(sig$gene_symbol, sig$gene_id)))),
  by_event_type = counts,
  top_events = top_list(sig, 10),
  top_by_event_type = setNames(lapply(types, function(t) top_list(sig[sig$event_type == t, ], 3)), types),
  heatmap = heat_made,
  splicing_vs_expression = integration,
  plots = list.files(out_dir, pattern = "\\.png$")
)
write_json_file(summary, file.path(out_dir, "splicing_summary.json"))
write_step_provenance(run_dir, STEP, started,
                      params_used = params[c("fdr_cutoff", "delta_psi_cutoff", "min_reads_per_sample", "n_label_events", "de_run",
                                             "n_heatmap_events", "plot_palette", "group_colours", "up_colour", "down_colour")],
                      packages = c("ggplot2", "ggrepel", "pheatmap"),
                      inputs = file.path(run_dir, "S01_validate_rmats", "events_validated.csv"), outputs = list.files(out_dir))
log_info("Step S02 done.")
