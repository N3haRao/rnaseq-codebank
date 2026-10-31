#!/usr/bin/env Rscript
# =============================================================================
# I04_isoform_report.R  --  methods paragraph for DEI and DUI
# =============================================================================
# Writes <run_dir>/I04_isoform_report/methods.md from the recorded parameters,
# results and package versions.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "isoforms/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "R", "gtf.R"))
source(file.path(REPO_ROOT, "isoforms", "R", "isoform_utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "I04_isoform_report"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "Isoform step I04: methods text."))
if (is_blank(opt$run_dir)) fail("--run_dir is required.")
params <- load_isoform_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
dei <- read_json_file(file.path(run_dir, "I02_differential_isoform_expression", "dei_summary.json"))
dui <- read_json_file(file.path(run_dir, "I03_differential_isoform_usage", "dui_summary.json"))
v <- function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "(not installed)"
ref <- params$reference_level; test <- params$test_level

methods <- c(
  "# Methods: transcript isoforms", "",
  paste0("Transcript abundances (", params$quant_type, ") for ", params$n_samples, " samples were imported with tximport ", v("tximport"),
         " (Soneson et al. 2015) at the transcript level and mapped to ", params$n_genes, " genes using the annotation ",
         "the transcriptome was built from."), "",
  paste0("**Differential isoform expression (DEI).** Transcript counts with effective-length offsets were analysed with DESeq2 ", v("DESeq2"),
         " (Love et al. 2014) using the design ", params$design_formula, ". Transcripts with at least ", params$min_count, " reads in at least ",
         params$min_samples, " samples were tested (Wald test, ", test, " vs ", ref, "; Benjamini-Hochberg adjustment); log2 fold changes were shrunk with ",
         dei$lfc_shrinkage_used, ". A transcript was called differentially expressed when ", dei$significance_rule, ": ",
         dei$n_significant, " transcripts (", dei$n_up, " higher and ", dei$n_down, " lower in ", test, ") in ",
         dei$n_genes_with_significant_transcripts, " genes."), "",
  if (!is.null(dui) && is.null(dui$status)) paste0(
    "**Differential isoform usage (DUI).** Changes in the share of each gene's expression made up by each transcript were tested following ",
    "Love, Soneson and Patro (2018): transcript counts scaled with countsFromAbundance = \"dtuScaledTPM\" were filtered (", dui$filter,
    ") and analysed with DEXSeq ", v("DEXSeq"), " (Anders et al. 2012) on the samples of the two groups, testing the condition-by-transcript ",
    "interaction within each gene. P-values were combined into gene-level q-values and controlled with stageR ", v("stageR"),
    " (Van den Berge et al. 2017): genes were screened and transcripts confirmed within significant genes at an overall FDR of ", params$dui_alpha,
    ". A transcript was called when ", dui$significance_rule, ": ", dui$n_transcripts_significant, " transcripts in ", dui$n_genes_significant,
    " of ", dui$n_genes_tested, " tested genes, including ", dui$n_isoform_switches, " genes with an isoform switch (one transcript gaining ",
    "and another losing share).")
  else paste0("Differential isoform usage was not tested: ", dui$reason %||% "no result", "."), "",
  paste0("Analysis code: RNA-seq code bank ", CODE_BANK_VERSION, ", isoform module; R ", R.version$major, ".", R.version$minor, "."), "",
  "References", "",
  "- Soneson C, Love MI, Robinson MD (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research 4:1521.",
  "- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15:550.",
  "- Anders S, Reyes A, Huber W (2012). Detecting differential usage of exons from RNA-seq data. Genome Research 22:2008-2017.",
  "- Van den Berge K, Soneson C, Robinson MD, Clement L (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151.",
  "- Love MI, Soneson C, Patro R (2018). Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification. F1000Research 7:952."
)
writeLines(methods, file.path(out_dir, "methods.md"))
write_step_provenance(run_dir, STEP, started, params_used = list(), outputs = list.files(out_dir))
log_info("Step I04 done: methods.md written.")
