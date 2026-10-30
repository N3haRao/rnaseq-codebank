#!/usr/bin/env Rscript
# =============================================================================
# U05_sashimi_plots.R  --  sashimi plots of the top significant splicing events
# =============================================================================
# R replacement for rMATS_ggsashimi.py + the plotting part of rmats_analysis.py.
# For each top event of a splicing-module run, it draws per group:
#   * mean read coverage across the event region (from the BAM files)
#   * splice-junction arcs, labelled with the mean number of reads per sample
#   * the annotated transcripts of the region (GTF), with the event's exon shaded
# Changes from the Python version: the plot window comes from the event's own
# exon coordinates (not from every numeric column of the table), nothing is
# plotted when no event is significant, and failures are reported.
#
# Reads BAM files with samtools (no Bioconductor packages needed); BAM files
# must be sorted and indexed (.bai), which STAR + samtools index give.
#
# Usage
#   Rscript later/scripts/U05_sashimi_plots.R --splicing_run runs/<splicing run> \
#       --comparison_dir results/rmats/WT_vs_KO --gtf ref/genes.gtf [--n_events 10]
#
# Outputs (in --out_dir, default <splicing_run>/sashimi/):
#   <rank>_<event_type>_<gene>.png, sashimi_summary.json
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "later/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "later", "R", "later_utils.R"))
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--splicing_run", type = "character", help = "Run folder of splicing/scripts/run_splicing.R (required)"),
  make_option("--comparison_dir", type = "character", help = "U04 comparison folder with b1.txt and b2.txt (required)"),
  make_option("--gtf", type = "character", help = "Annotation GTF (required)"),
  make_option("--n_events", type = "integer", default = 10L, help = "Top significant events to plot [default %default]"),
  make_option("--min_junction_reads", type = "double", default = 5, help = "Show junctions with at least this mean read count [default %default]"),
  make_option("--samtools", type = "character", help = "Path to samtools (default: on PATH)"),
  make_option("--out_dir", type = "character", help = "Output folder (default <splicing_run>/sashimi)")
), description = "Sashimi plots of significant splicing events."))
for (k in c("splicing_run", "comparison_dir", "gtf")) if (is_blank(opt[[k]])) fail("--", k, " is required.")
samtools <- need_tool("samtools", opt[["samtools"]])

res <- yaml::read_yaml(file.path(opt[["splicing_run"]], "S01_validate_rmats", "resolved_params.yaml"))
sig_path <- file.path(opt[["splicing_run"]], "S02_splicing_results", "splicing_events_significant.csv")
if (!file.exists(sig_path)) fail("No significant-events table in ", opt[["splicing_run"]], ". Run the splicing module first.")
sig <- data.table::fread(sig_path, data.table = FALSE, colClasses = list(character = c("gene_id", "gene_symbol", "chr")))
out_dir <- opt[["out_dir"]] %||% file.path(opt[["splicing_run"]], "sashimi")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!nrow(sig)) {
  write_json_file(list(n_plots = 0, note = "No significant events, so nothing was plotted."), file.path(out_dir, "sashimi_summary.json"))
  log_info("No significant events: nothing to plot.")
  quit(save = "no", status = 0)
}
top <- head(sig[order(sig$fdr, -abs(sig$delta_psi)), ], opt[["n_events"]])

bams <- list(strsplit(readLines(file.path(opt[["comparison_dir"]], "b1.txt"), n = 1), ",")[[1]],
             strsplit(readLines(file.path(opt[["comparison_dir"]], "b2.txt"), n = 1), ",")[[1]])
names(bams) <- c(res$group_1_name, res$group_2_name)
missing <- unlist(bams)[!file.exists(unlist(bams))]
if (length(missing)) fail("BAM file(s) not found: ", paste(missing, collapse = ", "))
for (b in unlist(bams)) if (!file.exists(paste0(b, ".bai")) && !file.exists(sub("\\.bam$", ".bai", b))) {
  if (run_tool(samtools, c("index", shQuote(b)), file.path(out_dir, "logs", "samtools_index.log")) != 0) {
    fail("Could not index ", b, " (is it sorted by coordinate, and is the folder writable?).")
  }
}
gtf <- read_gtf(opt[["gtf"]], features = "exon")

# ---- Reading one region from one BAM ----------------------------------------------------
coverage <- function(bam, chr, s, e) {
  out <- suppressWarnings(system2(samtools, c("depth", "-a", "-r", shQuote(sprintf("%s:%d-%d", chr, s, e)), shQuote(bam)),
                                  stdout = TRUE, stderr = FALSE))
  if (!length(out)) return(data.frame(pos = s:e, depth = 0))
  d <- data.table::fread(text = paste(out, collapse = "\n"), header = FALSE, data.table = FALSE)
  data.frame(pos = d$V2, depth = d$V3)
}
junctions <- function(bam, chr, s, e) {
  out <- suppressWarnings(system2(samtools, c("view", "-F", "0x904", shQuote(bam), shQuote(sprintf("%s:%d-%d", chr, s, e))),
                                  stdout = TRUE, stderr = FALSE))
  out <- out[grepl("^[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*N", out)]   # reads with a splice (N) in the CIGAR
  if (!length(out)) return(data.frame(start = integer(0), end = integer(0), n = integer(0)))
  f <- data.table::tstrsplit(out, "\t", fixed = TRUE, keep = c(4L, 6L))
  pos <- as.integer(f[[1]]); cig <- f[[2]]
  j <- do.call(rbind, lapply(seq_along(cig), function(i) {
    len <- as.integer(regmatches(cig[i], gregexpr("[0-9]+", cig[i]))[[1]])
    op <- regmatches(cig[i], gregexpr("[MIDNSHP=X]", cig[i]))[[1]]
    ref <- pos[i]; hits <- NULL
    for (k in seq_along(op)) {
      if (op[k] == "N") hits <- rbind(hits, c(ref - 1L, ref + len[k]))   # last base of donor exon, first base of acceptor exon
      if (op[k] %in% c("M", "D", "N", "=", "X")) ref <- ref + len[k]
    }
    hits
  }))
  j <- as.data.frame(j); names(j) <- c("start", "end")
  j <- j[j$start >= s & j$end <= e, ]
  if (!nrow(j)) return(data.frame(start = integer(0), end = integer(0), n = integer(0)))
  stats::aggregate(list(n = rep(1L, nrow(j))), by = list(start = j$start, end = j$end), FUN = sum)
}

# ---- One plot per event ------------------------------------------------------------------------
pal <- group_palette(c(res$group_2_name, res$group_1_name), list(plot_palette = "okabe_ito"))
made <- list()
for (r in seq_len(nrow(top))) {
  ev <- top[r, ]
  width <- ev$region_end - ev$region_start
  s <- max(1L, as.integer(ev$region_start - max(100, 0.1 * width)))
  e <- as.integer(ev$region_end + max(100, 0.1 * width))
  cov_df <- list(); jn_df <- list()
  for (g in names(bams)) {
    covs <- lapply(bams[[g]], coverage, chr = ev$chr, s = s, e = e)
    cm <- Reduce(function(a, b) merge(a, b, by = "pos", all = TRUE), lapply(seq_along(covs), function(i) setNames(covs[[i]], c("pos", paste0("d", i)))))
    cm[is.na(cm)] <- 0
    cov_df[[g]] <- data.frame(group = g, pos = cm$pos, depth = rowMeans(cm[, -1, drop = FALSE]))
    js <- do.call(rbind, lapply(bams[[g]], junctions, chr = ev$chr, s = s, e = e))
    if (!is.null(js) && nrow(js)) {
      js <- stats::aggregate(list(n = js$n), by = list(start = js$start, end = js$end), FUN = sum)
      js$n <- js$n / length(bams[[g]])
      js <- js[js$n >= opt[["min_junction_reads"]], ]
      if (nrow(js)) {
        dmap <- setNames(cov_df[[g]]$depth, cov_df[[g]]$pos)
        at <- function(x) { v <- unname(dmap[as.character(x)]); ifelse(is.na(v), 0, v) }
        js$y_start <- at(js$start)
        js$y_end <- at(js$end)
        js$group <- g
        js$up <- seq_len(nrow(js)) %% 2 == 1
        jn_df[[g]] <- js
      }
    }
  }
  cov_all <- do.call(rbind, cov_df)
  cov_all$group <- factor(cov_all$group, levels = names(bams))
  jn <- if (length(jn_df)) do.call(rbind, jn_df) else NULL
  ymax <- max(cov_all$depth, 1)
  gene <- display_name(ev$gene_symbol, ev$gene_id)
  p1 <- ggplot(cov_all, aes(pos, depth)) +
    annotate("rect", xmin = ev$event_start, xmax = ev$event_end, ymin = -Inf, ymax = Inf, fill = "#F0E442", alpha = 0.25) +
    geom_area(aes(fill = group), alpha = 0.85, show.legend = FALSE) +
    scale_fill_manual(values = pal) +
    facet_grid(group ~ ., switch = "y") +
    coord_cartesian(xlim = c(s, e), ylim = c(-0.7 * ymax, 1.6 * ymax)) +
    scale_y_continuous(breaks = pretty(c(0, ymax), 3)) +
    labs(title = sprintf("%s  %s  %s:%s-%s (%s)", gene, ev$event_type, ev$chr, format(ev$event_start, big.mark = ","),
                         format(ev$event_end, big.mark = ","), ev$strand),
         subtitle = sprintf("dPSI %+.2f (%s - %s), FDR %.2g. Mean coverage and junction reads per sample; shaded: alternative exon/region.",
                            ev$delta_psi, res$group_1_name, res$group_2_name, ev$fdr),
         x = NULL, y = "Mean coverage") +
    theme_codebank(10) + theme(strip.placement = "outside", axis.text.x = element_blank(), axis.ticks.x = element_blank())
  if (!is.null(jn) && nrow(jn)) {
    jn$group <- factor(jn$group, levels = names(bams))
    jn$lw <- 0.3 + 2.2 * jn$n / max(jn$n)
    for (u in c(TRUE, FALSE)) {
      d <- jn[jn$up == u, ]
      if (!nrow(d)) next
      p1 <- p1 + geom_curve(data = d, aes(x = start, xend = end, y = if (u) y_start else 0, yend = if (u) y_end else 0,
                                          colour = group, linewidth = lw),
                            curvature = if (u) -0.3 else 0.2, show.legend = FALSE)
      d$ylab <- if (u) pmax(d$y_start, d$y_end) + 0.35 * ymax else -0.45 * ymax
      p1 <- p1 + geom_label(data = d, aes(x = (start + end) / 2, y = ylab, label = round(n)), size = 2.6, label.padding = unit(0.1, "lines"))
    }
    p1 <- p1 + scale_colour_manual(values = pal) + scale_linewidth_identity()
  }
  ex <- gtf[gtf$chr == ev$chr & gtf$end >= s & gtf$start <= e, ]
  tx <- unique(ex$transcript_id)
  if (length(tx) > 8) tx <- tx[order(-table(ex$transcript_id)[tx])][1:8]
  ex <- ex[ex$transcript_id %in% tx, ]
  if (nrow(ex)) {
    ex$y <- match(ex$transcript_id, tx)
    ex$track <- "Transcripts"
    rng <- unique(ex[, c("transcript_id", "y", "track")])
    rng$start <- as.numeric(tapply(ex$start, ex$transcript_id, min)[rng$transcript_id])
    rng$end <- as.numeric(tapply(ex$end, ex$transcript_id, max)[rng$transcript_id])
    p2 <- ggplot() +
      annotate("rect", xmin = ev$event_start, xmax = ev$event_end, ymin = -Inf, ymax = Inf, fill = "#F0E442", alpha = 0.25) +
      geom_segment(data = rng, aes(x = pmax(start, s), xend = pmin(end, e), y = y, yend = y), colour = "grey40") +
      geom_rect(data = ex, aes(xmin = pmax(start, s), xmax = pmin(end, e), ymin = y - 0.3, ymax = y + 0.3), fill = "grey25") +
      scale_y_continuous(breaks = seq_along(tx), labels = tx) +
      facet_grid(track ~ ., switch = "y") +
      coord_cartesian(xlim = c(s, e)) + labs(x = sprintf("%s position", ev$chr), y = " ") +
      theme_codebank(9) + theme(strip.placement = "outside")
  } else {
    p2 <- ggplot() + annotate("text", x = (s + e) / 2, y = 0, label = "No annotated transcripts in this region") +
      coord_cartesian(xlim = c(s, e)) + theme_void()
  }
  file <- file.path(out_dir, sprintf("%02d_%s_%s.png", r, ev$event_type, gsub("[^A-Za-z0-9._-]", "_", gene)))
  h_cov <- 2.2 * length(bams); h_tx <- 0.6 + 0.3 * max(1, length(tx))
  png(file, width = 10, height = h_cov + h_tx, units = "in", res = 200)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = grid::unit(c(h_cov, h_tx), "null"))))
  g1 <- ggplotGrob(p1); g2 <- ggplotGrob(p2)
  if (length(g1$widths) == length(g2$widths)) {  # same layout: give both the same column widths so the x axes line up
    w <- grid::unit.pmax(g1$widths, g2$widths)
    g1$widths <- w; g2$widths <- w
  }
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1)); grid::grid.draw(g1); grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1)); grid::grid.draw(g2); grid::popViewport()
  invisible(dev.off())
  made[[length(made) + 1]] <- list(file = basename(file), event_id = ev$event_id, gene = gene, region = sprintf("%s:%d-%d", ev$chr, s, e),
                                   n_junctions_shown = if (is.null(jn)) 0L else nrow(jn))
  log_info("Saved ", basename(file))
}
write_json_file(list(n_plots = length(made), min_junction_reads = opt[["min_junction_reads"]], samtools = tool_version(samtools),
                     plots = made), file.path(out_dir, "sashimi_summary.json"))
log_info(sprintf("%d sashimi plot(s) -> %s", length(made), out_dir))
