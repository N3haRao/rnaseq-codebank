#!/usr/bin/env Rscript
# =============================================================================
# MJ04_majiq_deltapsi.R  --  change in PSI between the groups; voila tables
# =============================================================================
# majiq deltapsi -psi1 <reference>.psicov -psi2 <test>.psicov -n <reference> <test>
#   --changing-threshold <dpsi_threshold>   (dPSI = test - reference)
# voila tsv ... --show-all --show-read-counts            -> voila_all.tsv
# voila tsv ... --changing-between-group-dpsi <t> --probability-threshold <p>   -> voila_changing.tsv
# voila modulize ... sgc/ --changing-between-group-dpsi <t> --probability-changing-threshold <p>
#   --decomplexify-psi-threshold <d>      -> modulize/   (with run_modulize)
#
# Outputs (in <run_dir>/MJ04_majiq_deltapsi/): deltapsi.dpsicov, deltapsi.majiq.tsv,
# voila_all.tsv, voila_changing.tsv, modulize/, logs
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "majiq/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "majiq", "R", "majiq_utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "MJ04_majiq_deltapsi"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character"))))
params <- load_majiq_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
sg <- file.path(run_dir, "MJ02_majiq_build", "build", "splicegraph.zarr")
q <- file.path(run_dir, "MJ03_majiq_quantify")
ref <- params$ref_group; test <- params$test_group
lic <- params$majiq_license_file
dpsicov <- file.path(out_dir, "deltapsi.dpsicov")
code <- run_majiq_tool(params$majiq_bin_dir, "majiq",
                       c("deltapsi", "--license", lic, "-psi1", file.path(q, "psicov", paste0(ref, ".psicov")),
                         "-psi2", file.path(q, "psicov", paste0(test, ".psicov")), "-n", ref, test, "--splicegraph", sg,
                         "--output-tsv", file.path(out_dir, "deltapsi.majiq.tsv"), "--output-voila", dpsicov,
                         "--changing-threshold", params$dpsi_threshold, "--nthreads", params$threads,
                         "--logger", file.path(out_dir, "deltapsi_internal.log")),
                       file.path(out_dir, "majiq_deltapsi.log"))
if (code != 0) fail("majiq deltapsi failed (see majiq_deltapsi.log).")
code <- run_majiq_tool(params$majiq_bin_dir, "voila",
                       c("--license", lic, "tsv", sg, dpsicov, "-f", file.path(out_dir, "voila_all.tsv"),
                         "--show-all", "--show-read-counts", "-j", params$threads),
                       file.path(out_dir, "voila_tsv_all.log"))
if (code != 0) fail("voila tsv failed (see voila_tsv_all.log).")
code <- run_majiq_tool(params$majiq_bin_dir, "voila",
                       c("--license", lic, "tsv", sg, dpsicov, "-f", file.path(out_dir, "voila_changing.tsv"),
                         "--changing-between-group-dpsi", params$dpsi_threshold,
                         "--probability-threshold", params$probability_threshold, "-j", params$threads),
                       file.path(out_dir, "voila_tsv_changing.log"))
if (code != 0) fail("voila tsv (changing) failed (see voila_tsv_changing.log).")
mod_status <- "not run"
if (isTRUE(params$run_modulize)) {
  code <- run_majiq_tool(params$majiq_bin_dir, "voila",
                         c("--license", lic, "modulize", sg, dpsicov, file.path(q, "sgc"), "-d", file.path(out_dir, "modulize"),
                           "--overwrite", "--changing-between-group-dpsi", params$dpsi_threshold,
                           "--probability-changing-threshold", params$probability_threshold,
                           "--decomplexify-psi-threshold", params$decomplexify_psi_threshold, "-j", params$threads),
                         file.path(out_dir, "voila_modulize.log"))
  mod_status <- if (code == 0) "done" else "failed"
  if (code != 0) log_warn("voila modulize failed (see voila_modulize.log); the LSV tables are still complete.")
}
write_json_file(list(modulize = mod_status), file.path(out_dir, "deltapsi_status.json"))
write_step_provenance(run_dir, STEP, started, params_used = params[c("dpsi_threshold", "probability_threshold", "decomplexify_psi_threshold", "run_modulize")],
                      outputs = list.files(out_dir))
log_info("Step MJ04 done (modulize: ", mod_status, ").")
