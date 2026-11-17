#!/usr/bin/env Rscript
# =============================================================================
# SU03_suppa_diffsplice.R  --  test the change in PSI between the groups
# =============================================================================
# suppa.py diffSplice -m empirical (-gc: gene-level correction of p-values;
# -pa only for paired samples). The reference group goes first, so
# dPSI = PSI(test) - PSI(reference), the same sign as rMATS with --b1 = test.
# Significant: p <= pvalue_cutoff and |dPSI| >= dpsi_cutoff.
#
# Outputs (in <run_dir>/SU03_suppa_diffsplice/): diffsplice.dpsi / .psivec (SUPPA),
# suppa_events_all.csv, suppa_events_significant.csv, suppa_genes_significant.csv,
# suppa_summary.json, event_counts.png, dpsi_volcano.png, psi_scatter.png
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "suppa/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "suppa", "R", "suppa_utils.R"))
suppressPackageStartupMessages({ library(optparse); library(ggplot2) })

STEP <- "SU03_suppa_diffsplice"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character"))))
params <- load_suppa_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
d1 <- file.path(run_dir, "SU01_validate_suppa"); d2 <- file.path(run_dir, "SU02_suppa_events_psi")
cmd <- strsplit(params$suppa_command, " ", fixed = TRUE)[[1]]
ref <- params$reference_level; test <- params$test_level
set.seed(params$seed %||% 42)

# ---- diffSplice -----------------------------------------------------------------------------
if (isTRUE(params$paired) && params$group_sizes[[ref]] != params$group_sizes[[test]])
  fail("paired: true needs the same number of samples in both groups.")
psi <- file.path(d2, paste0(make.names(c(ref, test)), ".psi"))
tpm <- file.path(d1, c(params$tpm_files$reference, params$tpm_files$test))
prefix <- file.path(out_dir, "diffsplice")
args <- c("diffSplice", "-m", "empirical", "-i", file.path(d2, "events_all.ioe"), "-p", psi, "-e", tpm,
          "-a", params$area, "-l", params$lower_bound, "-th", params$tpm_threshold, "-o", prefix,
          if (isTRUE(params$gene_correction)) "-gc", if (isTRUE(params$paired)) "-pa")
code <- run_suppa(cmd, args, file.path(out_dir, "suppa_diffSplice.log"))
if (code != 0 || !file.exists(paste0(prefix, ".dpsi"))) fail("SUPPA diffSplice failed (see suppa_diffSplice.log).")

# ---- Tidy results ---------------------------------------------------------------------------
dp <- read_suppa_matrix(paste0(prefix, ".dpsi"))
if (ncol(dp) < 2) fail("Unexpected diffSplice output (", basename(prefix), ".dpsi).")
ev <- parse_event_ids(rownames(dp))
ev$dpsi <- dp[, 1]
ev$pvalue <- dp[, 2]
pv <- if (file.exists(paste0(prefix, ".psivec"))) read_suppa_matrix(paste0(prefix, ".psivec")) else NULL
meta <- data.table::fread(file.path(d1, "samples_validated.csv"), data.table = FALSE, colClasses = "character")
grp_mean <- function(g) {
  s <- intersect(meta$sample_id[meta[[params$condition]] == g], colnames(pv))
  if (is.null(pv) || !length(s)) return(rep(NA_real_, nrow(ev)))
  rowMeans(pv[match(ev$event_id, rownames(pv)), s, drop = FALSE], na.rm = TRUE)
}
ev$psi_reference <- grp_mean(ref)
ev$psi_test <- grp_mean(test)
tx <- data.table::fread(file.path(d1, "transcripts.csv"), data.table = FALSE)
gn <- tx[!duplicated(tx$gene_id), c("gene_id", "gene_name")]
ev$gene_name <- gn$gene_name[match(ev$gene_id, gn$gene_id)]
ev$gene_name[is.na(ev$gene_name) | ev$gene_name == ""] <- ev$gene_id[is.na(ev$gene_name) | ev$gene_name == ""]
ev$event_type_name <- unname(SUPPA_EVENT_NAMES[ev$event_type])
tested <- !is.na(ev$pvalue) & !is.na(ev$dpsi)
ev$significant <- tested & ev$pvalue <= params$pvalue_cutoff & abs(ev$dpsi) >= params$dpsi_cutoff
ev$direction <- ifelse(!ev$significant, "not significant", ifelse(ev$dpsi > 0, "more inclusion", "less inclusion"))
cols <- c("event_id", "gene_id", "gene_name", "event_type", "event_type_name", "chr", "strand", "event_start", "event_end",
          "coordinates", "psi_reference", "psi_test", "dpsi", "pvalue", "significant", "direction")
ev <- ev[, cols]
write_table(ev[order(ev$pvalue, -abs(ev$dpsi), na.last = TRUE), ], file.path(out_dir, "suppa_events_all.csv"))
sig <- ev[ev$significant, ]
sig <- sig[order(sig$pvalue, -abs(sig$dpsi)), ]
write_table(sig, file.path(out_dir, "suppa_events_significant.csv"))
genes <- if (nrow(sig)) {
  do.call(rbind, lapply(split(sig, sig$gene_id), function(d) data.frame(
    gene_id = d$gene_id[1], gene_name = d$gene_name[1], n_events = nrow(d),
    event_types = paste(sort(unique(d$event_type)), collapse = ";"),
    n_more_inclusion = sum(d$dpsi > 0), n_less_inclusion = sum(d$dpsi < 0),
    max_abs_dpsi = max(abs(d$dpsi)), min_pvalue = min(d$pvalue), stringsAsFactors = FALSE)))
} else data.frame(gene_id = character(0), gene_name = character(0), n_events = integer(0), event_types = character(0),
                  n_more_inclusion = integer(0), n_less_inclusion = integer(0), max_abs_dpsi = numeric(0), min_pvalue = numeric(0))
genes <- genes[order(genes$min_pvalue, -genes$max_abs_dpsi), , drop = FALSE]
write_table(genes, file.path(out_dir, "suppa_genes_significant.csv"))

# ---- Summary ---------------------------------------------------------------------------------
types <- intersect(SUPPA_EVENT_ORDER, unique(ev$event_type))
by_type <- lapply(types, function(t) list(event_type = t, name = SUPPA_EVENT_NAMES[[t]], tested = sum(tested & ev$event_type == t),
                                          more_inclusion = sum(sig$event_type == t & sig$dpsi > 0),
                                          less_inclusion = sum(sig$event_type == t & sig$dpsi < 0)))
rule <- sprintf("p <= %s%s and |dPSI| >= %s", params$pvalue_cutoff, if (isTRUE(params$gene_correction)) " (gene-corrected)" else "", params$dpsi_cutoff)
summ <- list(
  comparison = sprintf("%s vs %s", test, ref), method = "SUPPA2 diffSplice, empirical method",
  significance_rule = rule,
  sign_convention = sprintf("dPSI = PSI(%s) - PSI(%s); positive = more inclusion in %s", test, ref, test),
  n_events = nrow(ev), n_events_tested = sum(tested), n_significant = nrow(sig),
  n_up = sum(sig$dpsi > 0), n_down = sum(sig$dpsi < 0), n_genes_significant = nrow(genes),
  by_type = by_type,
  by_type_text = paste(vapply(by_type, function(b) sprintf("%s %d", b$event_type, b$more_inclusion + b$less_inclusion), ""), collapse = ", "),
  gene_correction = isTRUE(params$gene_correction), paired = isTRUE(params$paired),
  top_events = head(sig[, c("gene_name", "event_type", "dpsi", "pvalue", "event_id")], 10)
)
write_json_file(summ, file.path(out_dir, "suppa_summary.json"))
log_info(sprintf("%d of %d events significant (%s): %d more and %d less inclusion in %s, %d genes.",
                 nrow(sig), sum(tested), rule, summ$n_up, summ$n_down, test, nrow(genes)))

# ---- Figures ---------------------------------------------------------------------------------
cols_dir <- direction_colours(params)
dir_cols <- c(`more inclusion` = cols_dir[["up"]], `less inclusion` = cols_dir[["down"]], `not significant` = "grey70")
if (length(by_type)) {
  bt <- do.call(rbind, lapply(by_type, function(b) data.frame(event_type = b$event_type,
                                                             direction = c("more inclusion", "less inclusion"),
                                                             n = c(b$more_inclusion, b$less_inclusion))))
  bt$event_type <- factor(bt$event_type, levels = types)
  p <- ggplot(bt, aes(event_type, ifelse(direction == "more inclusion", n, -n), fill = direction)) +
    geom_col(width = 0.7) + geom_hline(yintercept = 0, colour = "grey40") +
    scale_fill_manual(values = dir_cols) + scale_y_continuous(labels = abs) +
    labs(x = NULL, y = "significant events", fill = NULL,
         title = sprintf("SUPPA2: significant events, %s vs %s", test, ref), subtitle = rule) + theme_codebank()
  save_plot(p, file.path(out_dir, "event_counts.png"), 7, 4.5)
}
vt <- ev[tested, ]
if (nrow(vt)) {
  # SUPPA's empirical p-values can be exactly 0: draw those just above the smallest non-zero p
  floor_p <- if (any(vt$pvalue > 0)) min(vt$pvalue[vt$pvalue > 0]) / 2 else 1e-6
  vt$neglog10p <- -log10(pmax(vt$pvalue, floor_p))
  lab <- head(vt[vt$significant, ][order(vt$pvalue[vt$significant]), ], params$n_label_events %||% 10)
  p <- ggplot(vt, aes(dpsi, neglog10p, colour = direction)) +
    geom_point(size = 0.9, alpha = 0.7) + scale_colour_manual(values = dir_cols) +
    geom_vline(xintercept = c(-1, 1) * params$dpsi_cutoff, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(params$pvalue_cutoff), linetype = "dashed", colour = "grey50") +
    labs(x = sprintf("dPSI (%s - %s)", test, ref), y = "-log10 p", colour = NULL,
         caption = if (any(vt$pvalue == 0)) "p = 0 drawn at half the smallest non-zero p" else NULL,
         title = "SUPPA2: change in inclusion per event") + theme_codebank()
  if (nrow(lab) && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = lab, aes(label = paste(gene_name, event_type)), size = 3, colour = "black", max.overlaps = 30)
  }
  save_plot(p, file.path(out_dir, "dpsi_volcano.png"), 7, 5.5)
  if (any(!is.na(vt$psi_reference))) {
    p <- ggplot(vt[order(vt$significant), ], aes(psi_reference, psi_test, colour = direction)) +
      geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed") +
      geom_point(size = 0.8, alpha = 0.7) + scale_colour_manual(values = dir_cols) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = sprintf("mean PSI, %s", ref), y = sprintf("mean PSI, %s", test), colour = NULL, title = "PSI per event in each group") +
      theme_codebank()
    save_plot(p, file.path(out_dir, "psi_scatter.png"), 6.5, 5.5)
  }
}
write_step_provenance(run_dir, STEP, started, packages = c("ggplot2"),
                      params_used = params[c("pvalue_cutoff", "dpsi_cutoff", "gene_correction", "paired", "area", "lower_bound", "tpm_threshold")],
                      outputs = list.files(out_dir))
log_info("Step SU03 done.")
