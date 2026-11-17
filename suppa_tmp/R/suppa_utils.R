# =============================================================================
# suppa/R/suppa_utils.R  --  helpers for the SUPPA2 module
# =============================================================================
# Source after R/utils.R and R/gtf.R.
# =============================================================================

SUPPA_ROOT <- file.path(REPO_ROOT, "suppa")

SUPPA_STEPS <- list(
  list(name = "SU01_validate_suppa",   critical = TRUE),
  list(name = "SU02_suppa_events_psi", critical = TRUE),
  list(name = "SU03_suppa_diffsplice", critical = TRUE),
  list(name = "SU04_suppa_vs_rmats",   critical = FALSE),
  list(name = "SU05_suppa_report",     critical = FALSE)
)

SUPPA_LIST_PARAMS <- c("exclude_samples", "event_types")

# SUPPA event codes -> the names used on the site and in rMATS
SUPPA_EVENT_NAMES <- c(SE = "skipped exon", A5 = "alternative 5' splice site", A3 = "alternative 3' splice site",
                       MX = "mutually exclusive exons", RI = "retained intron", AF = "alternative first exon",
                       AL = "alternative last exon")
SUPPA_EVENT_ORDER <- names(SUPPA_EVENT_NAMES)

load_suppa_defaults <- function() yaml::read_yaml(file.path(SUPPA_ROOT, "config", "default_params.yaml"))

load_suppa_params <- function(run_dir, use_resolved = TRUE) {
  params <- load_suppa_defaults()
  p <- file.path(run_dir, "params.yaml")
  if (file.exists(p)) params <- modifyList(params, yaml::read_yaml(p) %||% list())
  r <- file.path(run_dir, "SU01_validate_suppa", "resolved_params.yaml")
  if (use_resolved) {
    if (!file.exists(r)) fail("This step needs step SU01 (input check) to have passed first.")
    params <- modifyList(params, yaml::read_yaml(r) %||% list())
  }
  for (k in SUPPA_LIST_PARAMS) if (k %in% names(params)) params[[k]] <- as_list_param(params[[k]])
  params$run_dir <- normalizePath(run_dir)
  params
}

# ---- Finding SUPPA ----------------------------------------------------------------------------
# In order: the suppa_command parameter, $SUPPA_COMMAND, a conda environment made with
# suppa/install/install_suppa.sh ($COPILOT_HOME/envs/suppa or $SUPPA_ENV), suppa.py on the PATH.
# Returns the command as a character vector (program first), or NULL.
find_suppa <- function(params = list()) {
  split_cmd <- function(x) strsplit(trimws(x), "\\s+")[[1]]
  usable <- function(cmd) length(cmd) && (file.exists(cmd[1]) || nzchar(Sys.which(cmd[1])))
  for (c in c(params$suppa_command, Sys.getenv("SUPPA_COMMAND"))) {
    if (!is_blank(c)) { cmd <- split_cmd(c); if (usable(cmd)) return(cmd) }
  }
  envs <- c(Sys.getenv("SUPPA_ENV"), if (nzchar(Sys.getenv("COPILOT_HOME"))) file.path(Sys.getenv("COPILOT_HOME"), "envs", "suppa"))
  for (e in envs[nzchar(envs)]) {
    s <- file.path(e, "bin", "suppa.py")
    if (file.exists(s)) return(c(file.path(e, "bin", "python"), s)[file.exists(c(file.path(e, "bin", "python"), s))])
  }
  for (n in c("suppa.py", "suppa")) if (nzchar(Sys.which(n))) return(unname(Sys.which(n)))
  NULL
}

# Run a SUPPA sub-command; its output goes to log_file. Returns the exit code.
run_suppa <- function(cmd, args, log_file) {
  log_info("Running: ", paste(c(basename(cmd), args), collapse = " "))
  code <- suppressWarnings(system2(cmd[1], c(cmd[-1], args), stdout = log_file, stderr = log_file))
  if (!identical(as.integer(code), 0L)) {
    tail_lines <- if (file.exists(log_file)) tail(readLines(log_file, warn = FALSE), 15) else character(0)
    cat(paste0("    | ", tail_lines, "\n"), sep = "")
  }
  as.integer(code)
}

# ---- SUPPA files ----------------------------------------------------------------------------------
# .tpm / .psi files: the header has the sample names only; each row starts with the transcript / event ID.
write_suppa_matrix <- function(mat, path) {
  writeLines(paste(colnames(mat), collapse = "\t"), path)
  data.table::fwrite(data.table::data.table(id = rownames(mat), mat), path, sep = "\t", col.names = FALSE,
                     quote = FALSE, append = TRUE)
}

read_suppa_matrix <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 2) return(matrix(numeric(0), 0, 0))
  hdr <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  body <- data.table::fread(text = lines[-1], sep = "\t", header = FALSE, data.table = FALSE, na.strings = c("nan", "NaN", "NA", ""))
  ids <- as.character(body[[1]])
  m <- as.matrix(body[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  # a .dpsi file has 2 header fields and 3 columns per row; .psi/.tpm have one header field per sample
  colnames(m) <- if (length(hdr) == ncol(m)) hdr else if (length(hdr) == ncol(m) + 1) hdr[-1] else paste0("V", seq_len(ncol(m)))
  rownames(m) <- ids
  m
}

# Event IDs look like  ENSG...;SE:chr1:200-400:500-800:+  (gene;type:chr:coordinates...:strand)
parse_event_ids <- function(ids) {
  gene <- sub(";.*$", "", ids)
  rest <- sub("^[^;]*;", "", ids)
  parts <- strsplit(rest, ":", fixed = TRUE)
  type <- vapply(parts, `[`, "", 1)
  chr <- vapply(parts, `[`, "", 2)
  strand <- vapply(parts, function(p) p[length(p)], "")
  coords <- vapply(parts, function(p) paste(p[3:(length(p) - 1)], collapse = ":"), "")
  nums <- lapply(strsplit(coords, "[:-]"), function(x) suppressWarnings(as.numeric(x)))
  data.frame(event_id = ids, gene_id = gene, event_type = type, chr = chr, strand = strand, coordinates = coords,
             event_start = vapply(nums, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE), 0),
             event_end = vapply(nums, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE), 0),
             stringsAsFactors = FALSE)
}

# Splice junctions of a SUPPA event ("a-b" pairs: last base of the upstream exon .. first base of the
# downstream exon, 1-based) as one sorted key; RI is matched on its retained intron only.
suppa_junction_key <- function(type, coordinates) {
  vapply(seq_along(type), function(i) {
    j <- regmatches(coordinates[i], gregexpr("[0-9]+-[0-9]+", coordinates[i]))[[1]]
    if (type[i] == "RI") j <- j[1]
    paste(sort(j), collapse = ";")
  }, "")
}
