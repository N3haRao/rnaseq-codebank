#!/usr/bin/env Rscript
# =============================================================================
# 02_qc.R  --  per-sample quality control
# =============================================================================
# What it does
#   Summarises each sample's sequencing depth and count distribution, and
#   flags samples that look unusual compared with the rest.
#
# Metrics per sample
#   library_size        total reads assigned to genes
#   detected_genes      genes with at least one read
#   pct_zero_genes      % of genes with zero reads
#   filtered_genes_pct  % of the sample's reads that fall in genes kept by the low-count filter
#
# Flags (warnings, never automatic removal)
#   low_library_size    fewer reads than min_library_size
#   small_library       fewer than 10% of the median library size
#   few_detected_genes  fewer than 50% of the median number of detected genes
#
# Outputs (in <run_dir>/02_qc/)
#   qc_metrics.csv          one row per sample
#   qc_summary.json         totals, filter result and flags (read by the agent)
#   library_sizes.png       reads per sample
#   detected_genes.png      detected genes per sample
#   logcpm_boxplot.png      distribution of log2(CPM + 1) per sample (genes passing the filter)
#   logcpm_density.png      the same as density curves
#
# Usage
#   Rscript scripts/02_qc.R --run_dir runs/my_run
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
})

STEP <- "02_qc"
started <- Sys.time()

option_list <- list(
  make_option("--run_dir", type = "character", help = "Run folder (required)"),
  make_option("--config", type = "character", help = "Parameter file (default: <run_dir>/params.yaml)"),
  make_option("--set", type = "character", help = "Extra overrides, e.g. \"min_library_size=5e6\""),
  make_option("--min_library_size", type = "double", help = "Flag samples with fewer reads"),
  make_option("--min_count", type = "integer", help = "Low-count filter: minimum reads"),
  make_option("--min_samples", type = "integer", help = "Low-count filter: minimum samples")
)
opt <- parse_args(OptionParser(option_list = option_list, description = "Step 02: sample QC."))
params <- load_params(opt)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
log_info("Step 02: QC")

inp <- read_validated_inputs(run_dir)
counts <- inp$counts
meta <- inp$meta
grp <- condition_factor(meta, params)
pal <- group_palette(levels(grp), params)

# ---- Metrics ----------------------------------------------------------------
keep <- low_count_filter(counts, params$min_count, params$min_samples)
lib <- colSums(counts)
detected <- colSums(counts > 0)
metrics <- data.frame(
  sample_id = colnames(counts),
  group = as.character(grp),
  library_size = as.numeric(lib),
  library_size_millions = round(as.numeric(lib) / 1e6, 2),
  detected_genes = as.integer(detected),
  pct_zero_genes = round(100 * colMeans(counts == 0), 1),
  filtered_genes_pct = round(100 * colSums(counts[keep, , drop = FALSE]) / lib, 1),
  stringsAsFactors = FALSE
)

med_lib <- median(lib); med_det <- median(detected)
flag_rows <- list()
add_flag <- function(sample, flag, message) {
  flag_rows[[length(flag_rows) + 1]] <<- list(sample_id = sample, flag = flag, message = message)
}
for (s in colnames(counts)) {
  if (lib[[s]] < params$min_library_size) {
    add_flag(s, "low_library_size",
             sprintf("%s has %.1f million reads (below the %.1f million threshold).",
                     s, lib[[s]] / 1e6, params$min_library_size / 1e6))
  }
  if (lib[[s]] < 0.1 * med_lib) {
    add_flag(s, "small_library",
             sprintf("%s has %.0f%% of the median library size; it may have failed or been under-sequenced.",
                     s, 100 * lib[[s]] / med_lib))
  }
  if (detected[[s]] < 0.5 * med_det) {
    add_flag(s, "few_detected_genes",
             sprintf("%s detects %d genes versus a median of %d.", s, detected[[s]], as.integer(med_det)))
  }
}
metrics$flags <- vapply(metrics$sample_id, function(s) {
  f <- vapply(Filter(function(x) x$sample_id == s, flag_rows), `[[`, character(1), "flag")
  paste(f, collapse = ";")
}, character(1))
write_table(metrics, file.path(out_dir, "qc_metrics.csv"))
for (fr in flag_rows) log_warn(fr$message)

# ---- Plots ------------------------------------------------------------------
metrics$sample_id <- factor(metrics$sample_id, levels = colnames(counts)[order(grp)])
metrics$group <- factor(metrics$group, levels = levels(grp))
plot_width <- max(6, 0.35 * ncol(counts) + 3)

p_lib <- ggplot(metrics, aes(sample_id, library_size_millions, fill = group)) +
  geom_col() +
  geom_hline(yintercept = params$min_library_size / 1e6, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = pal, name = params$condition) +
  labs(title = "Library size per sample", x = NULL, y = "Reads assigned to genes (millions)",
       caption = sprintf("Dashed line: %.1f million-read QC threshold", params$min_library_size / 1e6)) +
  theme_codebank() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_lib, file.path(out_dir, "library_sizes.png"), width = plot_width)

p_det <- ggplot(metrics, aes(sample_id, detected_genes, fill = group)) +
  geom_col() +
  scale_fill_manual(values = pal, name = params$condition) +
  labs(title = "Detected genes per sample", x = NULL, y = "Genes with at least 1 read") +
  theme_codebank() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_det, file.path(out_dir, "detected_genes.png"), width = plot_width)

cpm <- t(t(counts[keep, , drop = FALSE]) / lib * 1e6)
logcpm <- log2(cpm + 1)
long <- data.frame(
  sample_id = factor(rep(colnames(logcpm), each = nrow(logcpm)), levels = levels(metrics$sample_id)),
  group = factor(rep(as.character(grp), each = nrow(logcpm)), levels = levels(grp)),
  value = as.vector(logcpm)
)
p_box <- ggplot(long, aes(sample_id, value, fill = group)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
  scale_fill_manual(values = pal, name = params$condition) +
  labs(title = "Expression distribution per sample",
       subtitle = sprintf("%d genes passing the low-count filter", sum(keep)),
       x = NULL, y = "log2(CPM + 1)") +
  theme_codebank() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_box, file.path(out_dir, "logcpm_boxplot.png"), width = plot_width)

p_den <- ggplot(long, aes(value, group = sample_id, colour = group)) +
  geom_density(linewidth = 0.4) +
  scale_colour_manual(values = pal, name = params$condition) +
  labs(title = "Expression density per sample", x = "log2(CPM + 1)", y = "Density",
       subtitle = "Curves of the same shape mean the samples are comparable") +
  theme_codebank()
save_plot(p_den, file.path(out_dir, "logcpm_density.png"))

# ---- Summary ----------------------------------------------------------------
summary <- list(
  n_samples = ncol(counts),
  n_genes_total = nrow(counts),
  n_genes_detected_any = sum(rowSums(counts) > 0),
  low_count_filter = sprintf(">= %d reads in >= %d samples", params$min_count, params$min_samples),
  n_genes_pass_filter = sum(keep),
  library_size = list(min = min(lib), median = med_lib, max = max(lib)),
  detected_genes = list(min = min(detected), median = med_det, max = max(detected)),
  n_flags = length(flag_rows),
  flags = flag_rows,
  plots = list.files(out_dir, pattern = "\\.png$")
)
write_json_file(summary, file.path(out_dir, "qc_summary.json"))
log_info(sprintf("%d of %d genes pass the low-count filter; %d sample flag(s).",
                 sum(keep), nrow(counts), length(flag_rows)))

write_step_provenance(run_dir, STEP, started,
                      params_used = params[c("min_count", "min_samples", "min_library_size")],
                      packages = c("ggplot2"),
                      inputs = file.path(run_dir, "01_validate_inputs", c("counts_validated.csv", "metadata_validated.csv")),
                      outputs = list.files(out_dir))
log_info("Step 02 done.")
