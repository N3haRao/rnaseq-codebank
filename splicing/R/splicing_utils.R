# =============================================================================
# splicing/R/splicing_utils.R  --  helpers for the alternative-splicing module
# =============================================================================
# Source after R/utils.R. Knows the layout of rMATS output tables
# (<EVENT>.MATS.JCEC.txt / <EVENT>.MATS.JC.txt) for the five event types.
# =============================================================================

SPLICING_ROOT <- file.path(REPO_ROOT, "splicing")

SPLICING_STEPS <- list(
  list(name = "S01_validate_rmats",      critical = TRUE),
  list(name = "S02_splicing_results",    critical = TRUE),
  list(name = "S03_splicing_enrichment", critical = FALSE),
  list(name = "S04_splicing_report",     critical = FALSE)
)

EVENT_TYPES <- c("SE", "MXE", "A3SS", "A5SS", "RI")

EVENT_LABELS <- c(
  SE   = "Skipped exon",
  MXE  = "Mutually exclusive exons",
  A3SS = "Alternative 3' splice site",
  A5SS = "Alternative 5' splice site",
  RI   = "Retained intron"
)

# Coordinate columns of each event type (rMATS, 0-based starts). The first pair
# is the alternative part of the event (the exon or region whose use changes).
EVENT_COORDS <- list(
  SE   = c("exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"),
  MXE  = c("1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd",
           "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"),
  A3SS = c("longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"),
  A5SS = c("longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"),
  RI   = c("riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
)

RMATS_COMMON <- c("GeneID", "geneSymbol", "chr", "strand", "IJC_SAMPLE_1", "SJC_SAMPLE_1",
                  "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncFormLen", "SkipFormLen", "PValue", "FDR",
                  "IncLevel1", "IncLevel2", "IncLevelDifference")

load_splicing_defaults <- function() {
  yaml::read_yaml(file.path(SPLICING_ROOT, "config", "default_params.yaml"))
}

SPLICING_LIST_PARAMS <- c("event_types", "gene_set_collections", "gmt_files", "group_colours",
                          "group_1_samples", "group_2_samples")

# Parameters of a splicing run: defaults < <run_dir>/params.yaml < resolved values from S01
load_splicing_params <- function(run_dir, use_resolved = TRUE) {
  params <- load_splicing_defaults()
  p <- file.path(run_dir, "params.yaml")
  if (file.exists(p)) params <- modifyList(params, yaml::read_yaml(p) %||% list())
  r <- file.path(run_dir, "S01_validate_rmats", "resolved_params.yaml")
  if (use_resolved) {
    if (!file.exists(r)) fail("This step needs step S01 (rMATS input check) to have passed first.")
    params <- modifyList(params, yaml::read_yaml(r) %||% list())
  }
  for (k in SPLICING_LIST_PARAMS) if (k %in% names(params)) params[[k]] <- as_list_param(params[[k]])
  params$run_dir <- normalizePath(run_dir)
  params
}

# rMATS writes some columns twice (ID) and quotes gene names; read it as text and clean up.
read_rmats_table <- function(path) {
  df <- data.table::fread(path, data.table = FALSE, check.names = FALSE, colClasses = "character",
                          showProgress = FALSE, na.strings = c("NA", ""))
  names(df) <- make.unique(trimws(names(df)))
  for (k in intersect(c("GeneID", "geneSymbol"), names(df))) df[[k]] <- gsub('^"|"$', "", df[[k]])
  df
}

# "12,5,NA" -> c(12, 5, NA)
split_numbers <- function(x) {
  lapply(strsplit(ifelse(is.na(x), "", x), ",", fixed = TRUE),
         function(v) suppressWarnings(as.numeric(v)))
}

n_values <- function(x) vapply(split_numbers(x), length, integer(1))

mean_of_list <- function(x) vapply(split_numbers(x), function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE), numeric(1))

# One rMATS table -> the module's standard columns (one row per event)
standardise_events <- function(df, event_type) {
  coords <- EVENT_COORDS[[event_type]]
  cm <- sapply(coords, function(k) suppressWarnings(as.numeric(df[[k]])))
  if (is.null(dim(cm))) cm <- matrix(cm, nrow = 1, dimnames = list(NULL, coords))
  # 0-based half-open starts in rMATS -> 1-based closed coordinates
  starts <- grepl("Start_0base$|ES$", coords)
  cm1 <- cm
  cm1[, starts] <- cm1[, starts] + 1
  id_col <- if ("ID" %in% names(df)) "ID" else names(df)[1]
  sym <- df$geneSymbol
  sym[is.na(sym) | sym %in% c("NA", "")] <- NA_character_
  ijc1 <- split_numbers(df$IJC_SAMPLE_1); sjc1 <- split_numbers(df$SJC_SAMPLE_1)
  ijc2 <- split_numbers(df$IJC_SAMPLE_2); sjc2 <- split_numbers(df$SJC_SAMPLE_2)
  reads_per_sample <- function(i, s) mapply(function(a, b) mean(a + b, na.rm = TRUE), i, s)
  out <- data.frame(
    event_type = event_type,
    event_id = paste0(event_type, "_", df[[id_col]]),
    gene_id = df$GeneID,
    gene_symbol = sym,
    chr = df$chr,
    strand = df$strand,
    event_start = cm1[, 1],
    event_end = cm[, 2],
    region_start = apply(cm1, 1, min, na.rm = TRUE),
    region_end = apply(cm, 1, max, na.rm = TRUE),
    mean_reads_1 = round(reads_per_sample(ijc1, sjc1), 2),
    mean_reads_2 = round(reads_per_sample(ijc2, sjc2), 2),
    psi_1 = round(mean_of_list(df$IncLevel1), 4),
    psi_2 = round(mean_of_list(df$IncLevel2), 4),
    delta_psi = suppressWarnings(as.numeric(df$IncLevelDifference)),
    pvalue = suppressWarnings(as.numeric(df$PValue)),
    fdr = suppressWarnings(as.numeric(df$FDR)),
    psi_values_1 = df$IncLevel1,
    psi_values_2 = df$IncLevel2,
    stringsAsFactors = FALSE
  )
  out$coordinates <- apply(cm1, 1, function(r) paste(r, collapse = ","))
  out
}

event_label <- function(ev) {
  sprintf("%s %s:%s-%s", display_name(ev$gene_symbol, ev$gene_id), ev$chr,
          format(ev$event_start, big.mark = ",", scientific = FALSE, trim = TRUE),
          format(ev$event_end, big.mark = ",", scientific = FALSE, trim = TRUE))
}

# Run-level step dispatch helpers (mirrors scripts/run_pipeline.R)
splicing_step_script <- function(name) file.path(SPLICING_ROOT, "scripts", paste0(name, ".R"))
