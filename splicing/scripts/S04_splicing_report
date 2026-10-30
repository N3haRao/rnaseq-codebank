#!/usr/bin/env Rscript
# =============================================================================
# S04_splicing_report.R  --  methods paragraph for the splicing analysis
# =============================================================================
# Writes <run_dir>/S04_splicing_report/methods.md from the recorded parameters,
# results and package versions (nothing is typed by hand).
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "splicing/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "splicing", "R", "splicing_utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "S04_splicing_report"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "Splicing step S04: methods text."))
if (is_blank(opt$run_dir)) fail("--run_dir is required.")
params <- load_splicing_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)

s <- read_json_file(file.path(run_dir, "S02_splicing_results", "splicing_summary.json"))
en <- read_json_file(file.path(run_dir, "S03_splicing_enrichment", "enrichment_summary.json"))
g1 <- params$group_1_name; g2 <- params$group_2_name
counts <- s$by_event_type
type_txt <- paste(sprintf("%d %s", counts$significant, counts$event_type), collapse = ", ")
rmats_version <- params$rmats_version %||% "(version not recorded; set rmats_version)"

methods <- c(
  "# Methods: alternative splicing",
  "",
  paste0(
    "Differential alternative splicing between ", g1, " (n = ", params$n_replicates_1, ") and ", g2,
    " (n = ", params$n_replicates_2, ") was analysed with rMATS ", rmats_version, " (Shen et al. 2014) for ",
    paste(params$event_types, collapse = ", "), " events, using ",
    if (params$count_type == "JCEC") "reads on splice junctions and on the alternative exon bodies (JCEC)" else "splice-junction reads only (JC)",
    ". rMATS estimates the percent spliced in (PSI) of each event per sample and tests the difference between groups ",
    "with a likelihood-ratio test; p-values were adjusted with the Benjamini-Hochberg method within each event type. ",
    "Events with fewer than ", params$min_reads_per_sample, " inclusion plus skipping reads per sample on average in either group were not ",
    "considered. An event was called significant when FDR <= ", params$fdr_cutoff, " and |dPSI| >= ", params$delta_psi_cutoff,
    " (dPSI = PSI in ", g1, " minus PSI in ", g2, "). Of ", s$n_events_tested, " events, ", s$n_events_passing_coverage,
    " had enough reads and ", s$n_significant, " were significant (", type_txt, "; ", s$n_up, " with more and ", s$n_down,
    " with less inclusion in ", g1, "), in ", s$n_genes_significant, " genes."),
  "",
  if (!is.null(en) && identical(en$status, "done")) paste0(
    "Pathway over-representation among genes with at least one significant event was tested with fgsea::fora ",
    "(fgsea ", en$package_versions$fgsea %||% "", "; hypergeometric test, Benjamini-Hochberg adjustment per collection) against all ",
    en$universe_size, " genes with an event that passed the read filter, using gene sets of ", params$gene_set_min_size, " to ",
    params$gene_set_max_size, " genes (", paste(vapply(en$collections, function(c) c$source, ""), collapse = "; "), ").")
  else paste0("Pathway enrichment was not run: ", (en$reason %||% "no result"), "."),
  "",
  paste0("Analysis code: RNA-seq code bank ", CODE_BANK_VERSION, ", splicing module; R ", R.version$major, ".", R.version$minor, "."),
  "",
  "References",
  "",
  "- Shen S, Park JW, Lu ZX, et al. (2014). rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. PNAS 111:E5593-E5601.",
  "- Korotkevich G, Sukhov V, Budin N, et al. (2021). Fast gene set enrichment analysis. bioRxiv 060012."
)
writeLines(methods, file.path(out_dir, "methods.md"))
write_step_provenance(run_dir, STEP, started, params_used = list(), outputs = list.files(out_dir))
log_info("Step S04 done: methods.md written.")
