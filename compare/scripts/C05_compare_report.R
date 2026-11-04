#!/usr/bin/env Rscript
# =============================================================================
# C05_compare_report.R  --  methods paragraph for the cross-study comparison
# =============================================================================
# Writes <run_dir>/C05_compare_report/methods.md from the recorded parameters,
# results and package versions (nothing is typed by hand).
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "compare/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "compare", "R", "compare_utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "C05_compare_report"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "Comparison step C05: methods text."))
if (is_blank(opt[["run_dir"]])) fail("--run_dir is required.")
params <- load_compare_params(opt[["run_dir"]])
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)

st <- read_json_file(file.path(run_dir, "C01_harmonise_studies", "studies.json"))
gc <- read_json_file(file.path(run_dir, "C02_compare_genes", "gene_comparison.json"))
rr <- read_json_file(file.path(run_dir, "C03_rrho", "rrho_summary.json"))
pw <- read_json_file(file.path(run_dir, "C04_compare_pathways", "pathway_comparison.json"))
s <- st$studies
n <- nrow(s)
study_txt <- paste(sprintf("%s (%s; %s, %s%s)", s$label, s$comparison, s$organism,
                           ifelse(s$source == "run", "analysed with the code bank's DESeq2 pipeline", "published results table"),
                           ifelse(s$flipped, "; fold changes sign-reversed to match the other studies", "")), collapse = "; ")
pk <- pkg_versions(c("babelgene", "fgsea", "msigdbr", "ggplot2"))

methods <- c(
  "# Methods: cross-study comparison",
  "",
  paste0(
    "The differential-expression results of ", n, " studies were compared: ", study_txt, ". ",
    "Each study was analysed separately, with its own normalisation and design; raw counts were not merged across studies. ",
    "Genes were matched by ", st$match_method,
    if (grepl("ortholog", st$match_method)) paste0(" (babelgene ", pk$babelgene %||% "", "; ortholog pairs supported by at least ",
                                                  params$min_ortholog_support, " databases)") else "",
    ". A gene was called significant in a study when ", st$significance_rule,
    ", using the same rule for every study. For each pair of studies, agreement of log2 fold changes was measured with ",
    "Pearson and Spearman correlation over the genes both studies tested, and the overlap of significant genes with the same ",
    "direction (up in both, down in both) was tested with one-sided Fisher's exact tests, using the genes tested in both studies as the universe."),
  "",
  if (!is.null(rr)) paste0(
    "Threshold-free agreement was assessed with the stratified rank-rank hypergeometric overlap test (RRHO2; Plaisier et al. 2010, ",
    "Cahill et al. 2018): genes were ranked by ", st$rank_by, " in each study and the overlap of the two rankings was tested ",
    "with the hypergeometric distribution at every pair of thresholds, separately for genes changing in the same and in opposite ",
    "directions; p values were ", if (identical(params$rrho_correction, "BY")) "adjusted with the Benjamini-Yekutieli method over the whole map" else "not adjusted", ".") else NULL,
  "",
  if (!is.null(pw) && identical(pw$status, "done")) paste0(
    "Pathways were compared using each study's pre-ranked GSEA results (fgsea ", pk$fgsea %||% "", ", ",
    paste(pw$collections_compared, collapse = ", "), " gene sets from MSigDB via msigdbr ", pk$msigdbr %||% "",
    "; ", params$gene_set_min_size, " to ", params$gene_set_max_size, " genes per set). Normalised enrichment scores were correlated ",
    "(Spearman) and pathways with padj < ", pw$padj_cutoff, " in both studies were counted by direction.") else
    if (!is.null(pw)) paste0("Pathways were not compared: ", pw$reason %||% "no result", ".") else NULL,
  "",
  paste0("Analysis code: RNA-seq code bank ", CODE_BANK_VERSION, ", compare module; R ", R.version$major, ".", R.version$minor, "."),
  "",
  "References",
  "",
  "- Plaisier SB, Taschereau R, Wong JA, Graeber TG. Rank-rank hypergeometric overlap: identification of statistically significant overlap between gene-expression signatures. Nucleic Acids Res. 2010;38:e169.",
  "- Cahill KM, Huo Z, Tseng GC, Logan RW, Seney ML. Improved identification of concordant and discordant gene expression signatures using an updated rank-rank hypergeometric overlap approach. Sci Rep. 2018;8:9588.",
  "- Korotkevich G, Sukhov V, Budin N, et al. Fast gene set enrichment analysis. bioRxiv. 2021. doi:10.1101/060012.",
  "- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15:550."
)
methods <- methods[!vapply(methods, is.null, logical(1))]
writeLines(unlist(methods), file.path(out_dir, "methods.md"))
log_info("Wrote methods.md")
write_step_provenance(run_dir, STEP, started, params_used = list(), packages = names(pk), outputs = list.files(out_dir))
