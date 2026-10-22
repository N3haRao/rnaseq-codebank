#!/usr/bin/env Rscript
# =============================================================================
# 07_report.R  --  methods paragraph and single-file HTML report
# =============================================================================
# What it does
#   1. Writes a methods paragraph built only from what this run recorded
#      (resolved parameters, step summaries and package versions), so the text
#      always matches the analysis that was actually run.
#   2. Renders templates/Report.Rmd into one self-contained HTML file with the
#      summary, warnings, all figures, top genes, pathways, methods and software.
#
#   The methods text does not depend on pandoc. If pandoc is missing, methods.md
#   is still written and the step finishes with a warning instead of failing.
#   Steps that did not run (e.g. enrichment skipped) are left out of the text.
#
# Outputs (in <run_dir>/07_report/)
#   methods.md    methods paragraph and references (Markdown)
#   report.html   full report (needs pandoc)
#
# Usage
#   Rscript scripts/07_report.R --run_dir runs/my_run
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "scripts/x"), ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "07_report"
started <- Sys.time()

option_list <- list(
  make_option("--run_dir", type = "character", help = "Run folder (required)"),
  make_option("--config", type = "character", help = "Parameter file (default: <run_dir>/params.yaml)"),
  make_option("--set", type = "character", help = "Extra overrides")
)
opt <- parse_args(OptionParser(option_list = option_list, description = "Step 07: methods text and HTML report."))
params <- load_params(opt)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
log_info("Step 07: report")

sdir <- function(step, file) file.path(run_dir, step, file)
val <- read_json_file(sdir("01_validate_inputs", "validation_report.json"))
qc  <- read_json_file(sdir("02_qc", "qc_summary.json"))
ex  <- read_json_file(sdir("03_exploratory", "exploratory_summary.json"))
de  <- read_json_file(sdir("04_differential_expression", "de_summary.json"))
en  <- read_json_file(sdir("06_enrichment", "enrichment_summary.json"))
if (is.null(de)) fail("Step 04 results not found (de_summary.json). The report needs steps 01-05.")

# Package versions as recorded by the steps that ran (not the versions loaded now)
pk <- list()
for (f in list.files(file.path(run_dir, "provenance"), pattern = "\\.json$", full.names = TRUE)) {
  p <- read_json_file(f)
  if (length(p$packages)) pk <- modifyList(pk, as.list(p$packages))
}
v <- function(pkg) if (is.null(pk[[pkg]])) "" else paste0(" version ", pk[[pkg]])

# Code-bank commit. provenance.json is written by run_pipeline.R after the last
# step, so it does not exist yet; ask git directly.
git_commit <- local({
  out <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(REPO_ROOT), "rev-parse", "--short", "HEAD"),
                                           stdout = TRUE, stderr = FALSE)),
                  error = function(e) character(0))
  if (length(out) && is.null(attr(out, "status"))) out[1] else ""
})

# -----------------------------------------------------------------------------
# 1. Methods paragraph
# -----------------------------------------------------------------------------
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
r_ver <- sub("^R version ([0-9.]+).*$", "\\1", R.version.string)
covs <- as_list_param(params$covariates)
n_samples <- val$summary$n_samples %||% length(de$size_factors)
excluded <- as_list_param(params$exclude_samples)

shrink_txt <- switch(de$lfc_shrinkage_used,
  apeglm = sprintf("Log2 fold changes were shrunken with apeglm%s (Zhu et al. 2019)", v("apeglm")),
  ashr   = sprintf("Log2 fold changes were shrunken with ashr%s (Stephens 2017)", v("ashr")),
  normal = "Log2 fold changes were shrunken with the DESeq2 'normal' prior",
  none   = "Log2 fold changes are reported without shrinkage",
  sprintf("Log2 fold changes were shrunken (%s)", de$lfc_shrinkage_used))

paras <- character(0)
paras <- c(paras, paste0(
  sprintf("Gene-level read counts for %d samples", n_samples),
  if (length(excluded)) sprintf(" (%s excluded before analysis)", paste(excluded, collapse = ", ")) else "",
  sprintf(" and %s genes were analysed in R %s with the RNA-seq code bank version %s",
          fmt_n(de$n_genes_input), r_ver, CODE_BANK_VERSION),
  if (nzchar(git_commit)) sprintf(" (commit %s)", git_commit) else "",
  ". ",
  sprintf("Genes with at least %d reads in at least %d samples were kept (%s of %s genes). ",
          params$min_count, params$min_samples, fmt_n(de$n_genes_tested), fmt_n(de$n_genes_input)),
  sprintf("Differential expression between %s and %s (metadata column '%s', reference level %s) was tested with DESeq2%s (Love et al. 2014) ",
          de$test_level, de$reference_level, de$condition, de$reference_level, v("DESeq2")),
  sprintf("using the design %s", de$design_formula),
  if (length(covs)) sprintf(", which adjusts for %s", paste(covs, collapse = " and ")) else "",
  if (length(de$numeric_covariates_scaled)) sprintf(" (numeric covariates %s were centred and scaled)",
                                                    paste(de$numeric_covariates_scaled, collapse = ", ")) else "",
  ". Size factors were estimated with the median-of-ratios method and dispersions with a ",
  de$dispersion_fit_type, " fit. ",
  "P-values from the Wald test were adjusted with the Benjamini-Hochberg method, with DESeq2's independent filtering. ",
  shrink_txt, ". ",
  sprintf("Genes were called differentially expressed when %s, giving %s genes (%s higher and %s lower in %s).",
          de$significance_rule, fmt_n(de$n_significant), fmt_n(de$n_up), fmt_n(de$n_down), de$test_level)
))

if (!is.null(ex)) {
  paras <- c(paras, sprintf(paste0(
    "For quality control and sample clustering, counts were transformed with the DESeq2 variance-stabilising ",
    "transformation (blind to the design). Principal component analysis used the %d most variable genes; ",
    "PC1 and PC2 explained %.1f%% and %.1f%% of the variance. Euclidean distances between samples were ",
    "clustered hierarchically."),
    ex$n_top_variable_genes, ex$percent_variance$PC1, ex$percent_variance$PC2))
}

annot_pkg <- orgdb_for(params$organism)
# Only when the annotation package was installed when step 04 ran (it records its version)
if (!is.null(annot_pkg) && !is.null(pk[[annot_pkg]]) && params$gene_id_type %in% c("ensembl", "entrez", "symbol")) {
  paras <- c(paras, sprintf("Genes were annotated with %s%s.", annot_pkg, v(annot_pkg)))
}

refs <- c(
  "Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550.",
  if (de$lfc_shrinkage_used == "apeglm")
    "Zhu A, Ibrahim JG, Love MI (2019). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 35:2084-2092.",
  if (de$lfc_shrinkage_used == "ashr")
    "Stephens M (2017). False discovery rates: a new deal. *Biostatistics* 18:275-294."
)

if (!is.null(en) && identical(en$status, "done") && length(en$collections)) {
  sources <- vapply(en$collections, function(cc) cc$source %||% "", character(1))
  uses_msig <- any(grepl("msigdb", sources, ignore.case = TRUE))
  paras <- c(paras, paste0(
    sprintf("Pathway enrichment used fgsea%s (Korotkevich et al. 2021) with the gene-set %s %s",
            v("fgsea"), if (length(sources) == 1) "collection" else "collections",
            paste(sprintf("%s (%s)", names(sources), sources), collapse = ", ")),
    if (uses_msig) sprintf("; MSigDB collections were accessed through msigdbr%s (Liberzon et al. 2015)", v("msigdbr")) else "",
    sprintf(". Gene sets with %s were tested. ", en$methods$gene_set_size),
    "Over-representation analysis was run separately on the up- and down-regulated genes (hypergeometric test, ",
    sprintf("background: the %s tested genes with a gene symbol), ", fmt_n(en$universe_size)),
    "and pre-ranked gene-set enrichment analysis was run on all tested genes ranked by the DESeq2 Wald statistic. ",
    sprintf("P-values were adjusted with the Benjamini-Hochberg method within each collection, and gene sets with adjusted p < %g were reported.",
            en$methods$padj_cutoff)
  ))
  refs <- c(refs, "Korotkevich G, Sukhov V, Budin N, et al. (2021). Fast gene set enrichment analysis. *bioRxiv* 060012.")
  if (uses_msig) refs <- c(refs, "Liberzon A, Birger C, Thorvaldsdottir H, et al. (2015). The Molecular Signatures Database hallmark gene set collection. *Cell Systems* 1:417-425.")
} else if (!is.null(en)) {
  paras <- c(paras, sprintf("Pathway enrichment was not run: %s", en$reason %||% "see the run log."))
}

methods_md <- c("# Methods", "", paste(paras, collapse = "\n\n"), "",
                "## References", "", paste0("- ", refs), "",
                sprintf("*Written automatically from run %s on %s.*", basename(run_dir), format(Sys.Date())))
methods_path <- file.path(out_dir, "methods.md")
writeLines(methods_md, methods_path)
log_info("Wrote methods.md")

# -----------------------------------------------------------------------------
# 2. HTML report
# -----------------------------------------------------------------------------
template <- file.path(REPO_ROOT, "templates", "Report.Rmd")
report_path <- file.path(out_dir, "report.html")
report_note <- "written"
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  report_note <- "skipped: the rmarkdown package is not installed"
} else if (!rmarkdown::pandoc_available()) {
  report_note <- "skipped: pandoc is not installed (methods.md and all results are still available)"
} else if (!file.exists(template)) {
  report_note <- paste("skipped: template not found:", template)
} else {
  # Render in a temporary folder so knitr's intermediate files never land in the run.
  work <- tempfile("report_")
  dir.create(work)
  file.copy(template, file.path(work, "Report.Rmd"))
  ok <- tryCatch({
    with_logged_warnings(rmarkdown::render(
      file.path(work, "Report.Rmd"), output_file = "report.html", output_dir = out_dir,
      intermediates_dir = work, knit_root_dir = work, envir = new.env(), quiet = TRUE,
      params = list(run_dir = run_dir, repo_root = REPO_ROOT, code_bank_commit = git_commit)))
    TRUE
  }, error = function(e) { log_warn("Rendering the HTML report failed: ", conditionMessage(e)); FALSE })
  unlink(work, recursive = TRUE)
  if (!ok) report_note <- "failed: see logs/07_report.log"
}
if (report_note == "written") log_info("Wrote report.html") else log_warn("HTML report ", report_note)

write_step_provenance(run_dir, STEP, started,
                      params_used = list(),
                      packages = c("rmarkdown", "knitr"),
                      inputs = c(sdir("04_differential_expression", "de_summary.json"),
                                 sdir("06_enrichment", "enrichment_summary.json")),
                      outputs = list.files(out_dir),
                      notes = list(report_html = report_note,
                                   pandoc = if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available())
                                     as.character(rmarkdown::pandoc_version()) else NULL))
# A failed render is not a failed step: the methods text and all results exist.
log_info("Step 07 done.")
