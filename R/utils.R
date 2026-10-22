# =============================================================================
# R/utils.R  --  shared helpers for every script in the code bank
# =============================================================================
# Each analysis script sources this file first. Nothing here performs an
# analysis on its own; it holds the plumbing that every step needs:
#
#   * logging and clean failure messages
#   * loading and merging parameters (defaults < config < resolved < CLI)
#   * reading and writing tables, JSON and YAML
#   * the run-folder layout and per-step provenance records
#   * gene ID / organism detection and gene annotation
#   * the low-count filter (shared so every step filters identically)
#   * plot styling and a colour-blind-safe palette
#
# Conventions
#   - Step scripts are run as:  Rscript scripts/<step>.R --run_dir <run folder>
#   - Each step writes into <run folder>/<step name>/ (e.g. runs/x/04_differential_expression/)
#   - Each step writes <run folder>/provenance/<step name>.json
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
  library(data.table)
})

# REPO_ROOT (the top folder of the repository) must be set by the calling
# script before it sources this file. Every script does this in its first lines.
if (!exists("REPO_ROOT")) stop("Set REPO_ROOT before sourcing R/utils.R")

CODE_BANK_VERSION <- "0.1.0"

# The pipeline steps, in order. The output folder of each step has the same
# name as its script (without ".R"). `critical = FALSE` means the pipeline keeps
# going if that step fails (the core results are still valid).
PIPELINE_STEPS <- list(
  list(name = "01_validate_inputs",        critical = TRUE),
  list(name = "02_qc",                     critical = TRUE),
  list(name = "03_exploratory",            critical = TRUE),
  list(name = "04_differential_expression", critical = TRUE),
  list(name = "05_de_plots",               critical = TRUE),
  list(name = "06_enrichment",             critical = FALSE),
  list(name = "07_report",                 critical = FALSE)
)

# Parameters that hold a list of values. On the command line they are given as
# comma-separated text ("batch,sex"); in YAML as a list ([batch, sex]).
LIST_PARAMS <- c("covariates", "gene_set_collections", "gmt_files", "exclude_samples")

# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

is_blank <- function(x) {
  is.null(x) || length(x) == 0 || (length(x) == 1 && (is.na(x) || identical(trimws(as.character(x)), "")))
}

timestamp_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------

log_msg <- function(..., level = "INFO") {
  cat(sprintf("[%s] %-5s %s\n", format(Sys.time(), "%H:%M:%S"), level, paste0(...)))
  flush.console()
}
log_info <- function(...) log_msg(..., level = "INFO")
log_warn <- function(...) log_msg(..., level = "WARN")

# Stop the step with a plain-language message and a non-zero exit code.
# The pipeline runner reads the exit code; the message ends up in the step log.
fail <- function(...) {
  log_msg(..., level = "ERROR")
  quit(save = "no", status = 1)
}

# Run an expression and log any R warnings it produces instead of dropping them.
with_logged_warnings <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    log_warn("R warning: ", conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}

# -----------------------------------------------------------------------------
# Paths and run-folder layout
# -----------------------------------------------------------------------------

# Resolve a path from a config file. Relative paths are tried against the
# current working directory first, then against the repository root.
resolve_path <- function(path) {
  if (is_blank(path)) return(NULL)
  path <- path.expand(path)
  if (file.exists(path)) return(normalizePath(path))
  alt <- file.path(REPO_ROOT, path)
  if (file.exists(alt)) return(normalizePath(alt))
  path  # return unchanged; the caller reports that it does not exist
}

step_dir <- function(run_dir, step) {
  d <- file.path(run_dir, step)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

# -----------------------------------------------------------------------------
# Reading and writing
# -----------------------------------------------------------------------------

# Read a table from .csv / .tsv / .txt (separator detected automatically) or
# .xlsx / .xls (first sheet). Column names are kept exactly as written.
read_table_file <- function(path) {
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path)))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      fail("Reading Excel files needs the 'readxl' package. Install it, or save the sheet as .csv.")
    }
    df <- as.data.frame(readxl::read_excel(path, sheet = 1), check.names = FALSE)
  } else {
    df <- data.table::fread(path, data.table = FALSE, check.names = FALSE,
                            strip.white = TRUE, showProgress = FALSE)
  }
  # fread names an empty first header "V1"; treat it as the ID column name.
  names(df) <- trimws(names(df))
  df
}

write_table <- function(df, path) {
  data.table::fwrite(df, path)
  invisible(path)
}

write_json_file <- function(x, path) {
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE,
                       null = "null", na = "null", digits = NA)
  invisible(path)
}

read_json_file <- function(path) {
  if (!file.exists(path)) return(NULL)
  jsonlite::read_json(path, simplifyVector = TRUE)
}

write_yaml_file <- function(x, path, header = NULL) {
  txt <- yaml::as.yaml(x)
  if (!is.null(header)) txt <- paste0(paste0("# ", header, collapse = "\n"), "\n", txt)
  writeLines(txt, path)
  invisible(path)
}

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

load_defaults <- function() {
  yaml::read_yaml(file.path(REPO_ROOT, "config", "default_params.yaml"))
}

# Turn "a,b" into c("a","b"); keep vectors/lists as character vectors.
as_list_param <- function(x) {
  if (is.null(x)) return(character(0))
  if (is.list(x)) x <- unlist(x)
  if (length(x) == 1 && is.character(x) && grepl(",", x)) x <- strsplit(x, ",")[[1]]
  x <- trimws(as.character(x))
  x[nzchar(x)]
}

normalise_params <- function(p) {
  for (k in LIST_PARAMS) if (k %in% names(p)) p[[k]] <- as_list_param(p[[k]])
  p
}

# Parse "--set 'a=1;b=[x, y]'" into a named list (values parsed as YAML so
# numbers, true/false and lists keep their types).
parse_set_string <- function(s) {
  if (is_blank(s)) return(list())
  out <- list()
  for (piece in strsplit(s, ";", fixed = TRUE)[[1]]) {
    piece <- trimws(piece)
    if (!nzchar(piece)) next
    kv <- regmatches(piece, regexpr("=", piece), invert = TRUE)[[1]]
    if (length(kv) != 2) fail("Could not read --set entry '", piece, "'. Use key=value.")
    out[[trimws(kv[1])]] <- yaml::yaml.load(trimws(kv[2]))
  }
  out
}

# Merge parameter sources for a step.
#   opt          : list from optparse (must contain run_dir; may contain config, set, ...)
#   use_resolved : read <run_dir>/01_validate_inputs/resolved_params.yaml (all steps after 01)
load_params <- function(opt, use_resolved = TRUE) {
  if (is_blank(opt$run_dir)) fail("--run_dir is required.")
  run_dir <- normalizePath(opt$run_dir, mustWork = FALSE)
  if (!dir.exists(run_dir)) fail("Run folder does not exist: ", run_dir)

  params <- load_defaults()

  config_path <- opt$config %||% file.path(run_dir, "params.yaml")
  if (file.exists(config_path)) {
    params <- modifyList(params, yaml::read_yaml(config_path) %||% list())
  } else if (!is.null(opt$config)) {
    fail("Config file not found: ", opt$config)
  }

  resolved_path <- file.path(run_dir, "01_validate_inputs", "resolved_params.yaml")
  if (use_resolved) {
    if (!file.exists(resolved_path)) {
      fail("This step needs step 01 (input validation) to have passed first. ",
           "Missing: ", resolved_path)
    }
    params <- modifyList(params, yaml::read_yaml(resolved_path) %||% list())
  }

  # Command-line overrides: every option that was given (non-NULL) except the
  # housekeeping ones.
  skip <- c("run_dir", "config", "help", "set")
  cli <- opt[setdiff(names(opt), skip)]
  cli <- cli[!vapply(cli, is.null, logical(1))]
  if (length(cli)) params <- modifyList(params, cli)
  params <- modifyList(params, parse_set_string(opt$set))

  params <- normalise_params(params)
  params$run_dir <- run_dir
  params
}

# -----------------------------------------------------------------------------
# Validated inputs (written by step 01, read by every later step)
# -----------------------------------------------------------------------------

read_validated_inputs <- function(run_dir) {
  d <- file.path(run_dir, "01_validate_inputs")
  cpath <- file.path(d, "counts_validated.csv")
  mpath <- file.path(d, "metadata_validated.csv")
  if (!file.exists(cpath) || !file.exists(mpath)) {
    fail("Validated inputs not found in ", d, ". Run step 01 first.")
  }
  cdf <- data.table::fread(cpath, data.table = FALSE, check.names = FALSE)
  counts <- as.matrix(cdf[, -1, drop = FALSE])
  rownames(counts) <- as.character(cdf[[1]])
  storage.mode(counts) <- "integer"

  meta <- data.table::fread(mpath, data.table = FALSE, check.names = FALSE,
                            colClasses = "character")
  rownames(meta) <- meta$sample_id
  # restore numeric columns that were numeric originally
  num_info <- read_json_file(file.path(d, "column_types.json"))
  if (!is.null(num_info)) {
    for (col in names(num_info)) {
      if (identical(num_info[[col]], "numeric") && col %in% names(meta)) {
        meta[[col]] <- as.numeric(meta[[col]])
      }
    }
  }
  stopifnot(identical(colnames(counts), rownames(meta)))
  list(counts = counts, meta = meta)
}

# Factor with the reference level first, test level second, others after.
condition_factor <- function(meta, params) {
  lv <- unique(as.character(meta[[params$condition]]))
  lv <- c(params$reference_level, params$test_level,
          setdiff(sort(lv), c(params$reference_level, params$test_level)))
  factor(as.character(meta[[params$condition]]), levels = lv)
}

# -----------------------------------------------------------------------------
# Low-count filter (identical everywhere)
# -----------------------------------------------------------------------------
# Keep genes with at least `min_count` reads in at least `min_samples` samples.
low_count_filter <- function(counts, min_count, min_samples) {
  rowSums(counts >= min_count) >= min_samples
}

# -----------------------------------------------------------------------------
# Gene IDs, organism, annotation
# -----------------------------------------------------------------------------

# Returns list(id_type = "ensembl" | "entrez" | "symbol" | "other", organism_guess = ...)
# Entrez IDs are plain numbers and do not show the organism: organism_guess is
# "unknown" and guess_organism_from_entrez() decides by looking the IDs up.
detect_gene_ids <- function(ids) {
  ids <- ids[!is.na(ids)]
  n <- length(ids)
  frac <- function(pattern) if (n == 0) 0 else mean(grepl(pattern, ids))
  if (frac("^ENSG[0-9]+") > 0.5)    return(list(id_type = "ensembl", organism_guess = "human"))
  if (frac("^ENSMUSG[0-9]+") > 0.5) return(list(id_type = "ensembl", organism_guess = "mouse"))
  if (frac("^FBgn[0-9]+") > 0.5)    return(list(id_type = "other",   organism_guess = "other"))
  if (frac("^ENS[A-Z]*G[0-9]+") > 0.5) return(list(id_type = "ensembl", organism_guess = "other"))
  # NCBI Entrez Gene IDs (e.g. GEO's NCBI-generated raw count tables)
  if (frac("^[0-9]+$") > 0.5)          return(list(id_type = "entrez",  organism_guess = "unknown"))
  # Gene symbols: human symbols are upper case (TP53), mouse are capitalised (Trp53)
  looks_symbol <- frac("^[A-Za-z][A-Za-z0-9.-]*$") > 0.5
  if (looks_symbol) {
    upper <- mean(ids == toupper(ids))
    title <- mean(grepl("^[A-Z][a-z0-9]", ids))
    if (upper > 0.6) return(list(id_type = "symbol", organism_guess = "human"))
    if (title > 0.6) return(list(id_type = "symbol", organism_guess = "mouse"))
  }
  list(id_type = "other", organism_guess = "other")
}

# Pick human or mouse for Entrez IDs by how many IDs each annotation package knows.
# Returns list(organism = "human" | "mouse" | "other", fraction = share of IDs matched,
#              checked = organisms whose annotation package was installed).
guess_organism_from_entrez <- function(ids) {
  best <- list(organism = "other", fraction = 0, checked = character(0))
  for (org in c("human", "mouse")) {
    pkg <- orgdb_for(org)
    if (!requireNamespace(pkg, quietly = TRUE)) next
    best$checked <- c(best$checked, org)
    db <- getExportedValue(pkg, pkg)
    f <- mean(ids %in% AnnotationDbi::keys(db, keytype = "ENTREZID"))
    if (f > best$fraction) { best$organism <- org; best$fraction <- f }
  }
  if (best$fraction < 0.5) best$organism <- "other"
  best
}

strip_ensembl_version <- function(ids) sub("^(ENS[A-Z]*G[0-9]+)\\.[0-9]+$", "\\1", ids)

orgdb_for <- function(organism) {
  switch(organism, human = "org.Hs.eg.db", mouse = "org.Mm.eg.db", NULL)
}

# Map gene IDs to symbol, full name and Entrez ID with the Bioconductor OrgDb
# package for the organism. Returns a data.frame with one row per input ID.
annotate_genes <- function(ids, organism, id_type) {
  out <- data.frame(gene_id = ids, gene_symbol = NA_character_,
                    gene_name = NA_character_, entrez_id = NA_character_,
                    stringsAsFactors = FALSE)
  pkg <- orgdb_for(organism)
  if (id_type == "entrez") out$entrez_id <- ids
  if (is.null(pkg) || !(id_type %in% c("ensembl", "entrez", "symbol"))) {
    if (id_type == "symbol") out$gene_symbol <- ids
    return(out)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    log_warn("Annotation package ", pkg, " is not installed; gene names will be missing. ",
             "Install it with BiocManager::install('", pkg, "').")
    if (id_type == "symbol") out$gene_symbol <- ids
    return(out)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  db <- get(pkg)
  keytype <- switch(id_type, ensembl = "ENSEMBL", entrez = "ENTREZID", symbol = "SYMBOL")
  map1 <- function(column) {
    suppressMessages(unname(AnnotationDbi::mapIds(db, keys = ids, keytype = keytype,
                                                  column = column, multiVals = "first")))
  }
  out$gene_symbol <- if (id_type == "symbol") ids else map1("SYMBOL")
  out$gene_name   <- map1("GENENAME")
  out$entrez_id   <- if (id_type == "entrez") ids else map1("ENTREZID")
  out
}

# A label for plots: symbol when available, otherwise the ID.
display_name <- function(symbol, id) ifelse(is.na(symbol) | symbol == "", id, symbol)

# -----------------------------------------------------------------------------
# Provenance
# -----------------------------------------------------------------------------

pkg_versions <- function(pkgs) {
  v <- vapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  }, character(1))
  as.list(v[!is.na(v)])
}

md5_of <- function(paths) {
  paths <- paths[file.exists(paths)]
  if (!length(paths)) return(list())
  as.list(setNames(unname(tools::md5sum(paths)), basename(paths)))
}

# Write <run_dir>/provenance/<step>.json describing exactly what the step did.
write_step_provenance <- function(run_dir, step, started, params_used,
                                  packages = character(0), inputs = character(0),
                                  outputs = character(0), notes = list()) {
  pdir <- file.path(run_dir, "provenance")
  dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
  rec <- list(
    step = step,
    code_bank_version = CODE_BANK_VERSION,
    started = format(started, "%Y-%m-%d %H:%M:%S %Z"),
    finished = timestamp_now(),
    seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
    r_version = R.version.string,
    parameters = params_used,
    packages = pkg_versions(unique(c("data.table", "jsonlite", "yaml", packages))),
    input_md5 = md5_of(inputs),
    outputs = as.list(basename(outputs)),
    notes = notes
  )
  write_json_file(rec, file.path(pdir, paste0(step, ".json")))
}

# -----------------------------------------------------------------------------
# Plot style
# -----------------------------------------------------------------------------

# Okabe-Ito colour-blind-safe palette (reference level gets grey-blue, test gets orange)
GROUP_COLOURS <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
                   "#56B4E9", "#F0E442", "#999999", "#000000")

group_palette <- function(levels) {
  cols <- rep(GROUP_COLOURS, length.out = length(levels))
  setNames(cols, levels)
}

theme_codebank <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "right")
}

save_plot <- function(plot, path, width = 7, height = 5) {
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = 300, bg = "white")
  log_info("Saved plot: ", basename(path))
  invisible(path)
}
