#!/usr/bin/env Rscript
# =============================================================================
# 03_exploratory.R  --  PCA, sample distances and outlier checks
# =============================================================================
# What it does
#   Applies DESeq2's variance-stabilising transformation (VST, blind = TRUE, so
#   the transform ignores the group labels, as recommended for QC), then:
#     * PCA on the most variable genes, coloured by group (and by each covariate)
#     * a heatmap of sample-to-sample distances
#     * checks each sample's position against the groups (possible mislabels and
#       outliers); these are flags for a person to review, never automatic removals
#
# Outlier logic (in PCA space, using the PCs that together explain >= 80% of the
# variance, between 2 and 5 PCs)
#   clusters_with_other_group : the sample is closer to the centre of another
#                               group than to the centre of its own group
#                               (own centre computed without the sample itself).
#                               Can mean a swapped label, or a strong batch effect.
#   far_from_own_group        : the sample's distance to its own group's centre
#                               is extreme (robust z > 3.5); only with >= 6 samples.
#
# Outputs (in <run_dir>/03_exploratory/)
#   pca_plot.png, pca_by_<covariate>.png, pca_scree.png
#   sample_distance_heatmap.png
#   pca_coordinates.csv     PC scores per sample (with metadata)
#   pca_variance.csv        % variance explained per PC
#   outlier_flags.csv       one row per sample with the checks above
#   exploratory_summary.json
#   vst_counts.csv.gz       VST matrix (genes passing the filter)
#   vsd.rds                 the DESeq2 transformed object (used by step 05)
#   size_factors.csv        DESeq2 normalisation factors per sample
#
# Usage
#   Rscript scripts/03_exploratory.R --run_dir runs/my_run [--n_top_variable_genes 1000]
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
  library(DESeq2)
  library(pheatmap)
})

STEP <- "03_exploratory"
started <- Sys.time()

option_list <- list(
  make_option("--run_dir", type = "character", help = "Run folder (required)"),
  make_option("--config", type = "character", help = "Parameter file (default: <run_dir>/params.yaml)"),
  make_option("--set", type = "character", help = "Extra overrides"),
  make_option("--n_top_variable_genes", type = "integer", help = "Genes used for PCA (default 500)")
)
opt <- parse_args(OptionParser(option_list = option_list, description = "Step 03: PCA and sample distances."))
params <- load_params(opt)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
log_info("Step 03: exploratory analysis")

inp <- read_validated_inputs(run_dir)
counts <- inp$counts
meta <- inp$meta
grp <- condition_factor(meta, params)
pal <- group_palette(levels(grp), params)
covs <- params$covariates

# ---- Variance-stabilising transformation -------------------------------------
keep <- low_count_filter(counts, params$min_count, params$min_samples)
dds <- DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE],
                              colData = data.frame(row.names = colnames(counts), group = grp),
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
vsd <- tryCatch(
  with_logged_warnings(vst(dds, blind = TRUE)),
  error = function(e) {
    log_warn("vst() failed (", conditionMessage(e), "); using varianceStabilizingTransformation() on all genes.")
    varianceStabilizingTransformation(dds, blind = TRUE)
  })
mat <- assay(vsd)
saveRDS(vsd, file.path(out_dir, "vsd.rds"))
vst_out <- data.frame(gene_id = rownames(mat), round(mat, 4), check.names = FALSE)
data.table::fwrite(vst_out, file.path(out_dir, "vst_counts.csv.gz"))
write_table(data.frame(sample_id = colnames(dds), size_factor = round(sizeFactors(dds), 4)),
            file.path(out_dir, "size_factors.csv"))

# ---- PCA ---------------------------------------------------------------------
n_top <- min(params$n_top_variable_genes, nrow(mat))
rv <- apply(mat, 1, var)
top <- order(rv, decreasing = TRUE)[seq_len(n_top)]
pca <- prcomp(t(mat[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
var_pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)
n_pc <- min(10, ncol(pca$x))
scores <- as.data.frame(pca$x[, seq_len(n_pc), drop = FALSE])
coords <- data.frame(meta, scores, check.names = FALSE)
coords$group <- grp
write_table(coords[, c("sample_id", setdiff(names(meta), "sample_id"), colnames(scores))],
            file.path(out_dir, "pca_coordinates.csv"))
write_table(data.frame(PC = paste0("PC", seq_along(var_pct)),
                       percent_variance = round(var_pct, 2),
                       cumulative_percent = round(cumsum(var_pct), 2)),
            file.path(out_dir, "pca_variance.csv"))

axis_lab <- function(i) sprintf("PC%d (%.1f%% of variance)", i, var_pct[i])
label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 30, show.legend = FALSE)
} else {
  geom_text(aes(label = sample_id), size = 3, vjust = -1, show.legend = FALSE)
}

# Plot options: what sets colour and shape (only the picture changes, not the PCA)
plot_column <- function(name, what) {
  if (is_blank(name) || identical(name, "none")) return(NULL)
  if (!name %in% names(meta)) {
    log_warn(what, " column '", name, "' is not in the sample table; ignored.")
    return(NULL)
  }
  name
}
colour_var <- plot_column(params$pca_colour_by, "pca_colour_by")
if (identical(colour_var, params$condition)) colour_var <- NULL
if (identical(params$pca_shape_by, "none")) {
  shape_var <- NULL
} else if (!is_blank(params$pca_shape_by)) {
  shape_var <- plot_column(params$pca_shape_by, "pca_shape_by")
  if (!is.null(shape_var) && length(unique(meta[[shape_var]])) > 6) {
    log_warn("pca_shape_by '", shape_var, "' has more than 6 values; shapes would be unreadable, so it is ignored.")
    shape_var <- NULL
  }
} else {
  shape_var <- NULL
  for (cv in covs) {
    if (!cv %in% params$numeric_covariates && length(unique(meta[[cv]])) <= 6) { shape_var <- cv; break }
  }
}
if (!isTRUE(params$pca_label_samples %||% TRUE)) label_layer <- NULL
colour_scale <- if (is.null(colour_var)) {
  scale_colour_manual(values = pal, name = params$condition)
} else if (is.numeric(meta[[colour_var]])) {
  scale_colour_viridis_c(name = colour_var)
} else {
  lv <- sort(unique(as.character(meta[[colour_var]])))
  scale_colour_manual(values = group_palette(lv, params), name = colour_var)
}
coords$.colour <- if (is.null(colour_var)) coords$group else if (is.numeric(meta[[colour_var]])) meta[[colour_var]] else factor(meta[[colour_var]])
p_pca <- ggplot(coords, aes(PC1, PC2, colour = .colour)) +
  { if (is.null(shape_var)) geom_point(size = 3.5)
    else geom_point(aes(shape = factor(.data[[shape_var]])), size = 3.5) } +
  label_layer +
  colour_scale +
  { if (!is.null(shape_var)) labs(shape = shape_var) } +
  labs(title = "PCA of samples",
       subtitle = sprintf("Top %d most variable genes, variance-stabilised counts", n_top),
       x = axis_lab(1), y = axis_lab(2)) +
  theme_codebank()
save_plot(p_pca, file.path(out_dir, "pca_plot.png"), width = 7.5, height = 5.5)

for (cv in covs) {
  is_numeric_cov <- cv %in% params$numeric_covariates
  p <- ggplot(coords, aes(PC1, PC2, colour = if (is_numeric_cov) .data[[cv]] else factor(.data[[cv]]))) +
    geom_point(size = 3.5) + label_layer +
    labs(title = sprintf("PCA coloured by %s", cv), colour = cv, x = axis_lab(1), y = axis_lab(2),
         subtitle = "If samples group by this variable, it is worth adjusting for") +
    theme_codebank()
  save_plot(p, file.path(out_dir, sprintf("pca_by_%s.png", make.names(cv))), width = 7.5, height = 5.5)
}

scree <- data.frame(PC = factor(paste0("PC", seq_len(n_pc)), levels = paste0("PC", seq_len(n_pc))),
                    pct = var_pct[seq_len(n_pc)])
p_scree <- ggplot(scree, aes(PC, pct)) + geom_col(fill = "#0072B2") +
  labs(title = "Variance explained by each principal component", x = NULL, y = "% of variance") +
  theme_codebank()
save_plot(p_scree, file.path(out_dir, "pca_scree.png"), width = 6, height = 4)

# ---- Sample distance heatmap ------------------------------------------------
d <- dist(t(mat))
dm <- as.matrix(d)
ann <- data.frame(row.names = colnames(mat), check.names = FALSE)
ann[[params$condition]] <- grp
for (cv in covs) ann[[cv]] <- if (cv %in% params$numeric_covariates) meta[[cv]] else factor(meta[[cv]])
ann_colours <- setNames(list(pal), params$condition)
hm_size <- max(5, 0.3 * ncol(mat) + 3)
pheatmap(dm, clustering_distance_rows = d, clustering_distance_cols = d,
         annotation_col = ann, annotation_colors = ann_colours,
         color = colorRampPalette(c("#08306B", "#6BAED6", "#F7FBFF"))(100),
         main = "Sample-to-sample distances (VST)", border_color = NA,
         filename = file.path(out_dir, "sample_distance_heatmap.png"),
         width = hm_size + 1.5, height = hm_size)
log_info("Saved plot: sample_distance_heatmap.png")

# ---- Outlier / mislabel checks ----------------------------------------------
cum <- cumsum(var_pct)
k <- max(2, min(5, which(cum >= 80)[1], ncol(pca$x)))
k <- min(k, ncol(pca$x))
X <- pca$x[, seq_len(k), drop = FALSE]
groups <- as.character(grp)
centre <- function(rows) colMeans(X[rows, , drop = FALSE])
flags <- data.frame(sample_id = rownames(X), group = groups, nearest_group = NA_character_,
                    dist_to_own_group = NA_real_, dist_to_nearest_other = NA_real_,
                    flag = "", message = "", stringsAsFactors = FALSE)
for (i in seq_len(nrow(X))) {
  own <- groups[i]
  own_rows <- setdiff(which(groups == own), i)
  if (length(own_rows) == 0) next
  d_own <- sqrt(sum((X[i, ] - centre(own_rows))^2))
  others <- setdiff(unique(groups), own)
  d_oth <- vapply(others, function(g) sqrt(sum((X[i, ] - centre(which(groups == g)))^2)), numeric(1))
  flags$dist_to_own_group[i] <- d_own
  flags$dist_to_nearest_other[i] <- min(d_oth)
  flags$nearest_group[i] <- if (min(d_oth) < d_own) others[which.min(d_oth)] else own
  if (min(d_oth) < d_own) {
    flags$flag[i] <- "clusters_with_other_group"
    flags$message[i] <- sprintf(
      "%s (labelled '%s') sits closer to the '%s' samples than to its own group in the PCA. This can mean a sample mix-up or mislabel, or that another factor (e.g. batch) drives the clustering. Check the lab records before removing it.",
      flags$sample_id[i], own, others[which.min(d_oth)])
  }
}
if (nrow(X) >= 6) {
  dd <- flags$dist_to_own_group
  m <- median(dd, na.rm = TRUE); s <- mad(dd, na.rm = TRUE)
  if (is.finite(s) && s > 0) {
    z <- (dd - m) / s
    far <- which(z > 3.5 & flags$flag == "")
    flags$flag[far] <- "far_from_own_group"
    flags$message[far] <- sprintf("%s is unusually far from the rest of its group in the PCA (robust z = %.1f). It may be an outlier worth checking.",
                                  flags$sample_id[far], z[far])
  }
}
flags$dist_to_own_group <- round(flags$dist_to_own_group, 2)
flags$dist_to_nearest_other <- round(flags$dist_to_nearest_other, 2)
write_table(flags, file.path(out_dir, "outlier_flags.csv"))
flagged <- flags[flags$flag != "", ]
for (m in flagged$message) log_warn(m)

# How much of PC1 / PC2 is explained by the groups and by each covariate (R^2)
r2 <- function(pc, v) {
  if (length(unique(v)) < 2) return(NA_real_)
  round(summary(lm(pc ~ v))$r.squared, 3)
}
assoc <- lapply(c(params$condition, covs), function(v) {
  vals <- if (v == params$condition) grp else meta[[v]]
  if (!(v %in% params$numeric_covariates)) vals <- factor(vals)
  list(variable = v, r2_PC1 = r2(pca$x[, 1], vals), r2_PC2 = r2(pca$x[, 2], vals))
})

summary <- list(
  transformation = "DESeq2 vst(blind = TRUE)",
  n_genes_used = nrow(mat),
  n_top_variable_genes = n_top,
  percent_variance = list(PC1 = round(var_pct[1], 1), PC2 = round(var_pct[2], 1),
                          PC3 = if (length(var_pct) >= 3) round(var_pct[3], 1) else NULL),
  pcs_used_for_outlier_check = k,
  pc_association = assoc,
  n_flagged_samples = nrow(flagged),
  flagged_samples = if (nrow(flagged)) lapply(seq_len(nrow(flagged)), function(i)
    list(sample_id = flagged$sample_id[i], flag = flagged$flag[i], message = flagged$message[i])) else list(),
  plots = list.files(out_dir, pattern = "\\.png$"),
  pca_plot = list(colour_by = colour_var %||% params$condition, shape_by = shape_var,
                  palette = params$plot_palette %||% "okabe_ito", group_colours = as.list(pal))
)
write_json_file(summary, file.path(out_dir, "exploratory_summary.json"))
log_info(sprintf("PC1 %.1f%%, PC2 %.1f%% of variance; %d sample(s) flagged.",
                 var_pct[1], var_pct[2], nrow(flagged)))

write_step_provenance(run_dir, STEP, started,
                      params_used = params[intersect(c("n_top_variable_genes", "min_count", "min_samples", "plot_palette",
                                                       "group_colours", "pca_colour_by", "pca_shape_by",
                                                       "pca_label_samples"), names(params))],
                      packages = c("DESeq2", "ggplot2", "pheatmap", "ggrepel"),
                      inputs = file.path(run_dir, "01_validate_inputs", c("counts_validated.csv", "metadata_validated.csv")),
                      outputs = list.files(out_dir))
log_info("Step 03 done.")
