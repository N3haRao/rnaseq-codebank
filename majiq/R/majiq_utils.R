# =============================================================================
# majiq/R/majiq_utils.R  --  helpers for the MAJIQ module
# =============================================================================
# Source after R/utils.R.
# =============================================================================

MAJIQ_ROOT <- file.path(REPO_ROOT, "majiq")

MAJIQ_STEPS <- list(
  list(name = "MJ01_validate_majiq", critical = TRUE),
  list(name = "MJ02_majiq_build",    critical = TRUE),
  list(name = "MJ03_majiq_quantify", critical = TRUE),
  list(name = "MJ04_majiq_deltapsi", critical = TRUE),
  list(name = "MJ05_majiq_results",  critical = FALSE),
  list(name = "MJ06_majiq_report",   critical = FALSE)
)

MAJIQ_LIST_PARAMS <- c("exclude_samples")

load_majiq_defaults <- function() yaml::read_yaml(file.path(MAJIQ_ROOT, "config", "default_params.yaml"))

load_majiq_params <- function(run_dir, use_resolved = TRUE) {
  params <- load_majiq_defaults()
  p <- file.path(run_dir, "params.yaml")
  if (file.exists(p)) params <- modifyList(params, yaml::read_yaml(p) %||% list())
  r <- file.path(run_dir, "MJ01_validate_majiq", "resolved_params.yaml")
  if (use_resolved) {
    if (!file.exists(r)) fail("This step needs step MJ01 (input check) to have passed first.")
    params <- modifyList(params, yaml::read_yaml(r) %||% list())
  }
  for (k in MAJIQ_LIST_PARAMS) if (k %in% names(params)) params[[k]] <- as_list_param(params[[k]])
  params$run_dir <- normalizePath(run_dir)
  params
}

# ---- Finding MAJIQ and its licence --------------------------------------------------------------
# The folder holding majiq (and voila, majiq-v3): majiq_bin_dir, $MAJIQ_BIN_DIR,
# $COPILOT_HOME/envs/majiq/bin (majiq/install/install_majiq.sh), or the folder of majiq on the PATH.
find_majiq_bin <- function(params = list()) {
  home <- Sys.getenv("COPILOT_HOME")
  cands <- c(params$majiq_bin_dir, Sys.getenv("MAJIQ_BIN_DIR"), if (nzchar(home)) file.path(home, "envs", "majiq", "bin"))
  for (d in cands) if (!is_blank(d) && file.exists(file.path(d, "majiq"))) return(normalizePath(d))
  w <- Sys.which("majiq")
  if (nzchar(w)) return(dirname(unname(w)))
  NULL
}

find_majiq_license <- function(params = list()) {
  home <- Sys.getenv("COPILOT_HOME")
  cands <- c(params$majiq_license_file, Sys.getenv("MAJIQ_LICENSE_FILE"))
  for (f in cands) if (!is_blank(f) && file.exists(f)) return(normalizePath(f))
  for (d in c(if (nzchar(home)) home, path.expand("~"), getwd())) {
    hit <- list.files(d, pattern = "^majiq_license.*", full.names = TRUE)
    if (length(hit)) return(normalizePath(hit[1]))
  }
  NULL
}

# A MAJIQ / voila command; output to log_file. Returns the exit code.
run_majiq_tool <- function(bin_dir, tool, args, log_file) {
  exe <- file.path(bin_dir, tool)
  log_info("Running: ", paste(c(tool, args), collapse = " "))
  code <- suppressWarnings(system2(exe, args, stdout = log_file, stderr = log_file))
  if (!identical(as.integer(code), 0L)) {
    cat(paste0("    | ", tail(readLines(log_file, warn = FALSE), 15), "\n"), sep = "")
  }
  as.integer(code)
}

# BAM of a sample: <dir>/<sample>.bam, <dir>/<sample>.sorted.bam, <dir>/<sample>/Aligned.sortedByCoord.out.bam
# (STAR in later/), or <dir>/<sample>/<sample>.bam
find_sample_bam <- function(bam_dir, sample) {
  c1 <- file.path(bam_dir, c(paste0(sample, ".bam"), paste0(sample, ".sorted.bam"),
                             file.path(sample, "Aligned.sortedByCoord.out.bam"), file.path(sample, paste0(sample, ".bam"))))
  hit <- c1[file.exists(c1)]
  if (length(hit)) hit[1] else NA_character_
}

# ---- voila tables -------------------------------------------------------------------------------------
# voila writes a JSON comment block ("# {...}") before the header; read the table after it.
read_voila_tsv <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  if (!length(lines)) return(NULL)
  data.table::fread(text = lines, sep = "\t", data.table = FALSE, colClasses = "character", fill = TRUE)
}

# Which column holds what, by name (voila 2 and 3 name some columns differently)
voila_columns <- function(cols, ref, test) {
  pick <- function(patterns) { for (p in patterns) { h <- grep(p, cols, ignore.case = TRUE, value = TRUE); if (length(h)) return(h[1]) }; NA_character_ }
  list(lsv_id = pick(c("^lsv_id$", "lsv.?id")),
       gene_id = pick(c("^gene_id$", "gene.?id")),
       gene_name = pick(c("^gene_name$", "gene.?name", "^gene$")),
       dpsi = pick(c("^mean_dpsi_per_lsv_junction$", "mean_dpsi", "e\\(dpsi\\)", "^dpsi")),
       prob = pick(c("^probability_changing$", "p\\(\\|dpsi\\|\\s*>=", "probability_changing", "prob.*chang")),
       psi_ref = pick(c(paste0("^", ref, "_mean_psi$"), paste0("^", make.names(ref), "_mean_psi$"), paste0("^", ref, ".*psi"))),
       psi_test = pick(c(paste0("^", test, "_mean_psi$"), paste0("^", make.names(test), "_mean_psi$"), paste0("^", test, ".*psi"))),
       lsv_type = pick(c("^lsv_type$", "lsv.?type")),
       chr = pick(c("^seqid$", "^chr$", "chrom")),
       strand = pick(c("^strand$")),
       junctions = pick(c("^junctions_coords$", "junctions.?coord")))
}

# "0.1;-0.3;0.2" -> numbers
split_num <- function(x) lapply(strsplit(ifelse(is.na(x), "", x), "[;,]"), function(v) suppressWarnings(as.numeric(v)))
