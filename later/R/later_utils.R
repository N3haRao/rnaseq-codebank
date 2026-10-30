# =============================================================================
# later/R/later_utils.R  --  helpers for the FASTQ / BAM scripts in later/
# =============================================================================
# Source after R/utils.R. These scripts start from sequencing files (FASTQ,
# BAM), which the website does not take yet; they are kept here for the Globus
# upload integration. Every external tool's exit code is checked: a failure
# stops the script with exit code 1 and a plain message (never a silent skip).
# =============================================================================

LATER_ROOT <- file.path(REPO_ROOT, "later")
source(file.path(REPO_ROOT, "R", "gtf.R"))

# Find a command-line tool or stop with a clear message.
need_tool <- function(name, path = NULL) {
  exe <- if (!is_blank(path)) path else Sys.which(name)
  if (is_blank(exe) || !nzchar(exe) || (!file.exists(exe) && !nzchar(Sys.which(exe)))) {
    fail(sprintf("'%s' was not found. Install it (e.g. conda install -c bioconda %s) or give its path.", name, name))
  }
  unname(exe)
}

tool_version <- function(exe, args = "--version") {
  out <- tryCatch(suppressWarnings(system2(exe, args, stdout = TRUE, stderr = TRUE)), error = function(e) character(0))
  out <- out[nzchar(trimws(out))]
  if (length(out)) trimws(out[length(out)]) else NA_character_
}

# Run a tool, writing its output to `log_file`. Returns the exit code (0 = ok).
# `timeout_s` uses the `timeout` command when it exists (Linux; `gtimeout` on macOS with coreutils).
run_tool <- function(exe, args, log_file, timeout_s = NULL, dry_run = FALSE) {
  cmd_txt <- paste(c(shQuote(exe), args), collapse = " ")
  if (dry_run) {
    log_info("[dry run] ", cmd_txt)
    return(0L)
  }
  dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
  writeLines(c(paste("#", timestamp_now()), paste("#", cmd_txt)), log_file)
  to <- Sys.which("timeout"); if (!nzchar(to)) to <- Sys.which("gtimeout")
  if (!is.null(timeout_s) && nzchar(to)) {  # system2 quotes the command itself; arguments are quoted by the caller
    args <- c(as.character(timeout_s), shQuote(exe), args)
    exe <- to
  }
  code <- suppressWarnings(system2(exe, args, stdout = log_file, stderr = log_file))
  if (!identical(as.integer(code), 0L)) {
    log_msg(sprintf("Command failed (exit code %s): %s", code, cmd_txt), level = "ERROR")
    tail_lines <- tail(readLines(log_file, warn = FALSE), 10)
    cat(paste0("    | ", tail_lines, "\n"), sep = "")
  }
  as.integer(code)
}

# ---- Sample sheet: one row per FASTQ file -------------------------------------------
# Columns (names configurable): file_name, sample, read_orientation (R1 / R2; SE or
# R1 alone = single-end), plus any sample characteristics (condition, batch, ...).
read_sample_sheet <- function(path, sample_col = "sample", file_col = "file_name", orient_col = "read_orientation") {
  df <- read_table_file(path)
  miss <- setdiff(c(sample_col, file_col), names(df))
  if (length(miss)) fail(sprintf("The sample sheet is missing column(s): %s. Found: %s.", paste(miss, collapse = ", "), paste(names(df), collapse = ", ")))
  for (k in intersect(c(sample_col, file_col, orient_col), names(df))) df[[k]] <- trimws(as.character(df[[k]]))
  if (!orient_col %in% names(df)) df[[orient_col]] <- "SE"
  df
}

# Pair R1 / R2 per sample. Returns data.frame(sample, fastq_1, fastq_2 (NA = single-end), <other columns>)
pair_fastqs <- function(df, fastq_dir, sample_col = "sample", file_col = "file_name", orient_col = "read_orientation") {
  resolve_fq <- function(f) {
    if (file.exists(f)) return(normalizePath(f))
    p <- file.path(fastq_dir %||% ".", f)
    if (file.exists(p)) return(normalizePath(p))
    NA_character_
  }
  problems <- character(0)
  rows <- lapply(split(df, df[[sample_col]]), function(g) {
    o <- toupper(g[[orient_col]])
    r1 <- g[[file_col]][o %in% c("R1", "1", "SE", "")]
    r2 <- g[[file_col]][o %in% c("R2", "2")]
    s <- g[[sample_col]][1]
    if (length(r1) != 1) { problems <<- c(problems, sprintf("%s: %d R1/SE files (need exactly 1)", s, length(r1))); return(NULL) }
    if (length(r2) > 1) { problems <<- c(problems, sprintf("%s: %d R2 files (need 0 or 1)", s, length(r2))); return(NULL) }
    f1 <- resolve_fq(r1); f2 <- if (length(r2)) resolve_fq(r2) else NA_character_
    if (is.na(f1)) problems <<- c(problems, sprintf("%s: file not found: %s", s, r1))
    if (length(r2) && is.na(f2)) problems <<- c(problems, sprintf("%s: file not found: %s", s, r2))
    extra <- g[1, setdiff(names(g), c(sample_col, file_col, orient_col)), drop = FALSE]
    data.frame(sample = s, fastq_1 = f1, fastq_2 = f2, extra, check.names = FALSE, stringsAsFactors = FALSE)
  })
  if (length(problems)) fail("Sample sheet problems:\n  ", paste(problems, collapse = "\n  "))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Salmon's expected_format -> strandedness used by rMATS (and featureCounts / STAR)
SALMON_STRAND <- c(ISR = "reverse", SR = "reverse", ISF = "forward", SF = "forward", IU = "unstranded", U = "unstranded",
                   OSR = "reverse", OSF = "forward", OU = "unstranded", MSR = "reverse", MSF = "forward", MU = "unstranded")
RMATS_LIBTYPE <- c(unstranded = "fr-unstranded", forward = "fr-secondstrand", reverse = "fr-firststrand")
