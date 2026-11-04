# =============================================================================
# compare/R/compare_utils.R  --  helpers for the cross-study comparison module
# =============================================================================
# Source after R/utils.R. The module compares the results of differential-
# expression analyses (one per study); it never merges raw counts.
# =============================================================================

COMPARE_ROOT <- file.path(REPO_ROOT, "compare")

COMPARE_STEPS <- list(
  list(name = "C01_harmonise_studies", critical = TRUE),
  list(name = "C02_compare_genes",     critical = TRUE),
  list(name = "C03_rrho",              critical = FALSE),
  list(name = "C04_compare_pathways",  critical = FALSE),
  list(name = "C05_compare_report",    critical = FALSE)
)

COMPARE_LIST_PARAMS <- c("gene_set_collections", "gmt_files")

load_compare_defaults <- function() {
  yaml::read_yaml(file.path(COMPARE_ROOT, "config", "default_params.yaml"))
}

# modifyList() drops unnamed list entries, so `studies` (a list of studies) is replaced, not merged
merge_compare_params <- function(base, new) {
  new <- new %||% list()
  out <- modifyList(base, new[setdiff(names(new), "studies")])
  if (!is.null(new$studies)) out$studies <- new$studies
  out
}

compare_step_script <- function(name) file.path(COMPARE_ROOT, "scripts", paste0(name, ".R"))

# Parameters of a comparison run: defaults < <run_dir>/params.yaml < resolved values from C01
load_compare_params <- function(run_dir, use_resolved = TRUE) {
  params <- load_compare_defaults()
  p <- file.path(run_dir, "params.yaml")
  if (file.exists(p)) params <- merge_compare_params(params, yaml::read_yaml(p))
  r <- file.path(run_dir, "C01_harmonise_studies", "resolved_params.yaml")
  if (use_resolved) {
    if (!file.exists(r)) fail("This step needs step C01 (harmonise studies) to have passed first.")
    params <- merge_compare_params(params, yaml::read_yaml(r))
  }
  for (k in COMPARE_LIST_PARAMS) if (k %in% names(params)) params[[k]] <- as_list_param(params[[k]])
  params$run_dir <- normalizePath(run_dir)
  params
}

# A label that is safe in file names: "GSE60450 (luminal)" -> "GSE60450_luminal"
safe_label <- function(x) {
  x <- gsub("[^A-Za-z0-9.-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "study")
}

# Studies as harmonised by C01: list(info = <studies.json>, tables = list(<label> = data.frame))
read_harmonised <- function(run_dir) {
  d <- file.path(run_dir, "C01_harmonise_studies")
  info <- read_json_file(file.path(d, "studies.json"))
  if (is.null(info)) fail("C01_harmonise_studies/studies.json not found: run step C01 first.")
  tables <- list()
  for (i in seq_len(nrow(info$studies))) {
    s <- info$studies[i, ]
    tables[[s$label]] <- data.table::fread(file.path(d, "studies", paste0(s$file_label, ".csv")), data.table = FALSE,
                                           colClasses = list(character = c("key", "gene")))
  }
  list(info = info, tables = tables)
}

# All pairs of studies, in the order given: list(c(a, b), ...)
study_pairs <- function(labels) {
  if (length(labels) < 2) return(list())
  cmb <- utils::combn(seq_along(labels), 2)
  lapply(seq_len(ncol(cmb)), function(j) labels[cmb[, j]])
}

pair_folder <- function(a, b) paste0(safe_label(a), "_vs_", safe_label(b))

# Find a column by any of several common names (case-insensitive). Returns the name or NULL.
find_column <- function(cols, candidates) {
  low <- tolower(cols)
  for (c in tolower(candidates)) {
    hit <- which(low == c)
    if (length(hit)) return(cols[hit[1]])
  }
  NULL
}

# ---- RRHO: rank-rank hypergeometric overlap (RRHO2's four-quadrant version) -------
# a, b: named numeric vectors (same genes), larger = more up. Both lists are ordered
# from most up to most down. For each pair of thresholds (i, j) the overlap is taken
# in the quadrant the thresholds fall in:
#   up-up      first i of A  with first j of B   (both up)
#   down-down  rest of A     with rest of B      (both down)
#   up-down    first i of A  with rest of B      (up in A, down in B)
#   down-up    rest of A     with first j of B
# and its over-enrichment p value is the hypergeometric upper tail (Plaisier et al.
# 2010; Cahill et al. 2018). Counts come from 2-D cumulative sums of the binned ranks.
rrho_map <- function(a, b, step = NULL, correction = "BY") {
  genes <- intersect(names(a), names(b))
  a <- a[genes]; b <- b[genes]
  n <- length(genes)
  if (n < 50) stop("too few shared genes for RRHO (", n, ")")
  if (is.null(step) || is.na(step)) step <- max(1L, as.integer(floor(n / 100)))
  step <- as.integer(step)
  ra <- rank(-a, ties.method = "first")
  rb <- rank(-b, ties.method = "first")
  th <- seq(step, n - 1L, by = step)
  k <- length(th)
  # bins: (0, th1], (th1, th2], ..., (th_k, n]
  breaks <- c(0, th, n)
  ba <- findInterval(ra, breaks, left.open = TRUE)
  bb <- findInterval(rb, breaks, left.open = TRUE)
  m <- k + 1L
  tab <- matrix(tabulate((ba - 1L) * m + bb, nbins = m * m), m, m, byrow = TRUE)   # [bin in A, bin in B]
  cum <- t(apply(apply(tab, 2, cumsum), 1, cumsum))                                # 2-D cumulative sum
  le <- cum[seq_len(k), seq_len(k), drop = FALSE]  # genes with rank_A <= th_i and rank_B <= th_j
  ti <- matrix(th, k, k)              # rows: threshold in A
  tj <- matrix(th, k, k, byrow = TRUE)
  n_up_a <- sum(a > 0); n_up_b <- sum(b > 0)
  q_a_up <- ti <= max(n_up_a, 1)
  q_b_up <- tj <= max(n_up_b, 1)
  # overlap size, size of the A part, size of the B part, for the cell's quadrant
  ov <- ifelse(q_a_up & q_b_up, le,
        ifelse(!q_a_up & !q_b_up, n - ti - tj + le,
        ifelse(q_a_up & !q_b_up, ti - le, tj - le)))
  sa <- ifelse(q_a_up, ti, n - ti)
  sb <- ifelse(q_b_up, tj, n - tj)
  p <- stats::phyper(ov - 1, sa, n - sa, sb, lower.tail = FALSE)
  p[sa == 0 | sb == 0] <- 1
  quadrant <- ifelse(q_a_up & q_b_up, "up_up", ifelse(!q_a_up & !q_b_up, "down_down",
              ifelse(q_a_up & !q_b_up, "up_down", "down_up")))
  padj <- if (identical(correction, "BY")) matrix(stats::p.adjust(p, "BY"), k, k) else p
  score <- -log10(pmax(padj, 1e-300))
  grid <- data.frame(threshold_a = as.vector(ti), threshold_b = as.vector(tj), quadrant = as.vector(quadrant),
                     overlap = as.vector(ov), size_a = as.vector(sa), size_b = as.vector(sb),
                     expected = signif(as.vector(sa * sb / n), 4), pvalue = as.vector(p),
                     padj = as.vector(padj), neg_log10_padj = as.vector(score))
  best <- lapply(c("up_up", "down_down", "up_down", "down_up"), function(q) {
    g <- grid[grid$quadrant == q, ]
    if (!nrow(g)) return(NULL)
    r <- g[which.max(g$neg_log10_padj), ]
    in_a <- if (q %in% c("up_up", "up_down")) ra <= r$threshold_a else ra > r$threshold_a
    in_b <- if (q %in% c("up_up", "down_up")) rb <= r$threshold_b else rb > r$threshold_b
    list(quadrant = q, threshold_a = r$threshold_a, threshold_b = r$threshold_b, overlap = r$overlap,
         expected = r$expected, padj = signif(r$padj, 3), neg_log10_padj = round(r$neg_log10_padj, 2),
         genes = genes[in_a & in_b])
  })
  names(best) <- c("up_up", "down_down", "up_down", "down_up")
  list(grid = grid, best = best, n_genes = n, step = step, n_up_a = n_up_a, n_up_b = n_up_b)
}
