#!/usr/bin/env Rscript
# =============================================================================
# MJ06_majiq_report.R  --  methods paragraph for the MAJIQ run
# =============================================================================
# Output: <run_dir>/MJ06_majiq_report/methods.md
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "majiq/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "majiq", "R", "majiq_utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "MJ06_majiq_report"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character"))))
params <- load_majiq_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
s <- read_json_file(file.path(run_dir, "MJ05_majiq_results", "majiq_summary.json"))
ref <- params$reference_level; test <- params$test_level

methods <- c(
  "# Methods: alternative splicing with MAJIQ", "",
  paste0("Local splicing variations (LSVs) were quantified with MAJIQ ", params$majiq_version, " (Vaquero-Garcia et al. 2016, 2023). ",
         "A splice graph was built from the annotation (", basename(params$gff3_used), ") and the aligned reads of ",
         params$n_samples, " samples (majiq build, strandness ", params$strandness, "; a junction had to pass MAJIQ's read filters in ",
         "at least half of the samples of one group). PSI coverage was computed per group (majiq psi-coverage) and the change in ",
         "PSI between ", test, " (", params$group_sizes[[test]], " samples) and ", ref, " (", params$group_sizes[[ref]],
         " samples) estimated with majiq deltapsi. ", s$sign_convention %||% "", ". An LSV was called changing when ",
         s$significance_rule %||% "", if (!is.null(s)) paste0(": ", s$n_lsvs_changing, " of ", s$n_lsvs_quantified, " LSVs in ",
         s$n_genes_changing, " genes") else "", "."),
  "",
  if (identical(s$modulize, "done")) paste0("Changing LSVs were grouped into splicing modules and classified into event types with voila modulize ",
                                            "(decomplexify PSI threshold ", params$decomplexify_psi_threshold, "): ", s$modulize_text, "."),
  "",
  paste0("Analysis code: RNA-seq code bank ", CODE_BANK_VERSION, ", MAJIQ module; R ", R.version$major, ".", R.version$minor,
         ". MAJIQ is used under its academic licence."), "",
  "References", "",
  "- Vaquero-Garcia J, Barrera A, Gazzara MR, Gonzalez-Vallinas J, Lahens NF, Hogenesch JB, Lynch KW, Barash Y (2016). A new view of transcriptome complexity and regulation through the lens of local splicing variations. eLife 5:e11752.",
  "- Vaquero-Garcia J, Aicher JK, Jewell S, Gazzara MR, Radens CM, Jha A, Norton SS, Lahens NF, Grant GR, Barash Y (2023). RNA splicing analysis using heterogeneous and large RNA-seq datasets. Nature Communications 14:1230."
)
writeLines(methods[!is.na(methods)], file.path(out_dir, "methods.md"))
write_step_provenance(run_dir, STEP, started, params_used = list(), outputs = list.files(out_dir))
log_info("Step MJ06 done: methods.md written.")
