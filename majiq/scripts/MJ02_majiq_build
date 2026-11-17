#!/usr/bin/env Rscript
# =============================================================================
# MJ02_majiq_build.R  --  splice graph from the annotation and the reads
# =============================================================================
# majiq build --license <lic> <gff3> experiments.tsv build --strandness <s> --nthreads <n>
# (build group = the sample's group: a junction must pass MAJIQ's filters in at least
# half of the samples of one group, MAJIQ's default min-experiments 0.5)
#
# Outputs (in <run_dir>/MJ02_majiq_build/): build/splicegraph.zarr, build/<sample>.sj,
# majiq_build.log
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "majiq/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "majiq", "R", "majiq_utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "MJ02_majiq_build"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character"))))
params <- load_majiq_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
build <- file.path(out_dir, "build")
code <- run_majiq_tool(params$majiq_bin_dir, "majiq",
                       c("build", "--license", params$majiq_license_file, params$gff3_used,
                         file.path(run_dir, "MJ01_validate_majiq", "experiments.tsv"), build,
                         "--strandness", params$strandness, "--nthreads", params$threads,
                         "--logger", file.path(out_dir, "majiq_build_internal.log")),
                       file.path(out_dir, "majiq_build.log"))
if (code != 0) fail("majiq build failed (see majiq_build.log).")
sj <- list.files(build, pattern = "\\.sj$")
if (!dir.exists(file.path(build, "splicegraph.zarr")) || !length(sj)) fail("majiq build wrote no splicegraph.zarr / .sj files.")
write_step_provenance(run_dir, STEP, started, params_used = params[c("strandness", "threads")],
                      outputs = list.files(out_dir), notes = list(majiq_version = params$majiq_version))
log_info(sprintf("Step MJ02 done: splice graph and %d junction files.", length(sj)))
