# =============================================================================
# isoforms/R/isoform_utils.R  --  helpers for the isoform module
# =============================================================================
# Source after R/utils.R and R/gtf.R.
# =============================================================================

ISOFORM_ROOT <- file.path(REPO_ROOT, "isoforms")

ISOFORM_STEPS <- list(
  list(name = "I01_validate_isoforms",              critical = TRUE),
  list(name = "I02_differential_isoform_expression", critical = TRUE),
  list(name = "I03_differential_isoform_usage",      critical = FALSE),
  list(name = "I04_isoform_report",                  critical = FALSE)
)

ISOFORM_LIST_PARAMS <- c("covariates", "exclude_samples", "group_colours")

load_isoform_defaults <- function() yaml::read_yaml(file.path(ISOFORM_ROOT, "config", "default_params.yaml"))

load_isoform_params <- function(run_dir, use_resolved = TRUE) {
  params <- load_isoform_defaults()
  p <- file.path(run_dir, "params.yaml")
  if (file.exists(p)) params <- modifyList(params, yaml::read_yaml(p) %||% list())
  r <- file.path(run_dir, "I01_validate_isoforms", "resolved_params.yaml")
  if (use_resolved) {
    if (!file.exists(r)) fail("This step needs step I01 (input check) to have passed first.")
    params <- modifyList(params, yaml::read_yaml(r) %||% list())
  }
  for (k in ISOFORM_LIST_PARAMS) if (k %in% names(params)) params[[k]] <- as_list_param(params[[k]])
  params$run_dir <- normalizePath(run_dir)
  params
}

quant_file_name <- function(type) switch(type, salmon = "quant.sf", kallisto = "abundance.tsv",
                                         stop("quant_type must be salmon or kallisto"))

# Sample table, tximport objects and transcript annotation written by I01
read_isoform_inputs <- function(run_dir) {
  d <- file.path(run_dir, "I01_validate_isoforms")
  list(meta = data.table::fread(file.path(d, "samples_validated.csv"), data.table = FALSE, colClasses = "character"),
       txi = readRDS(file.path(d, "txi_counts.rds")),
       txi_dtu = readRDS(file.path(d, "txi_dtu_scaled.rds")),
       tx = data.table::fread(file.path(d, "transcripts.csv"), data.table = FALSE, colClasses = "character"))
}

# DESeq2 / DEXSeq need syntactic names: model data plus a map back to the original names
model_frame <- function(meta, params) {
  cd <- data.frame(row.names = meta$sample_id)
  maps <- list()
  for (v in c(params$covariates, params$condition)) {
    lv <- if (v == params$condition) c(params$reference_level, params$test_level,
                                       setdiff(sort(unique(meta[[v]])), c(params$reference_level, params$test_level)))
          else sort(unique(meta[[v]]))
    map <- setNames(make.names(lv, unique = TRUE), lv)
    maps[[v]] <- map
    cd[[make.names(v)]] <- factor(unname(map[meta[[v]]]), levels = unname(map))
  }
  list(col_data = cd, maps = maps)
}
