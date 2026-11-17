#!/usr/bin/env Rscript
# =============================================================================
# MJ03_majiq_quantify.R  --  PSI coverage per group (and splice-graph coverage)
# =============================================================================
# majiq psi-coverage --license <lic> build/splicegraph.zarr psicov/<group>.psicov <sample>.sj ...
# majiq-v3 sg-coverage --license <lic> build/splicegraph.zarr sgc/<group>.sgc <sample>.sj ...
#   (only with run_modulize: voila modulize needs the .sgc files)
#
# Outputs (in <run_dir>/MJ03_majiq_quantify/): psicov/<group>.psicov, sgc/<group>.sgc, logs
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "majiq/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "majiq", "R", "majiq_utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "MJ03_majiq_quantify"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character"))))
params <- load_majiq_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
build <- file.path(run_dir, "MJ02_majiq_build", "build")
dir.create(file.path(out_dir, "psicov"), showWarnings = FALSE)
if (isTRUE(params$run_modulize)) dir.create(file.path(out_dir, "sgc"), showWarnings = FALSE)
groups <- list(params$samples$reference, params$samples$test)
names(groups) <- c(params$ref_group, params$test_group)
for (g in names(groups)) {
  sj <- file.path(build, paste0(groups[[g]], ".sj"))
  if (!all(file.exists(sj))) fail("Missing junction files from majiq build: ", paste(basename(sj[!file.exists(sj)]), collapse = ", "))
  code <- run_majiq_tool(params$majiq_bin_dir, "majiq",
                         c("psi-coverage", "--license", params$majiq_license_file, file.path(build, "splicegraph.zarr"),
                           file.path(out_dir, "psicov", paste0(g, ".psicov")), sj, "--nthreads", params$threads,
                           "--logger", file.path(out_dir, paste0("psicov_", g, "_internal.log"))),
                         file.path(out_dir, paste0("psicov_", g, ".log")))
  if (code != 0) fail("majiq psi-coverage failed for ", g, ".")
  if (isTRUE(params$run_modulize)) {
    code <- run_majiq_tool(params$majiq_bin_dir, "majiq-v3",
                           c("sg-coverage", "--license", params$majiq_license_file, file.path(build, "splicegraph.zarr"),
                             file.path(out_dir, "sgc", paste0(g, ".sgc")), sj, "--nthreads", params$threads,
                             "--logger", file.path(out_dir, paste0("sgc_", g, "_internal.log"))),
                           file.path(out_dir, paste0("sgc_", g, ".log")))
    if (code != 0) fail("majiq-v3 sg-coverage failed for ", g, ".")
  }
}
write_step_provenance(run_dir, STEP, started, params_used = params[c("threads", "run_modulize")], outputs = list.files(out_dir, recursive = TRUE))
log_info("Step MJ03 done.")
