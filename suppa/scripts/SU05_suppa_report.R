#!/usr/bin/env Rscript
# =============================================================================
# SU05_suppa_report.R  --  methods paragraph for the SUPPA2 run
# =============================================================================
# Output: <run_dir>/SU05_suppa_report/methods.md
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "suppa/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "suppa", "R", "suppa_utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "SU05_suppa_report"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "SUPPA step SU05: methods text."))
params <- load_suppa_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
ev <- read_json_file(file.path(run_dir, "SU02_suppa_events_psi", "events_summary.json"))
s <- read_json_file(file.path(run_dir, "SU03_suppa_diffsplice", "suppa_summary.json"))
cmp <- read_json_file(file.path(run_dir, "SU04_suppa_vs_rmats", "comparison_summary.json"))
ref <- params$reference_level; test <- params$test_level
type_names <- c(SE = "skipped exons", SS = "alternative 5' and 3' splice sites", MX = "mutually exclusive exons",
                RI = "retained introns", FL = "alternative first and last exons")

methods <- c(
  "# Methods: alternative splicing with SUPPA2", "",
  paste0("Alternative-splicing events (", paste(type_names[params$event_types], collapse = ", "), ") were derived from the annotation ",
         "with SUPPA2 generateEvents (Trincado et al. 2018; strict boundaries), giving ", ev$n_events, " events. Transcript ",
         "abundances (Salmon TPM) of ", params$n_samples, " samples were converted to percent spliced in (PSI) per event and ",
         "sample with psiPerEvent. Differential splicing between ", test, " and ", ref, " (", params$group_sizes[[test]], " and ",
         params$group_sizes[[ref]], " samples) was tested with diffSplice using the empirical method",
         if (isTRUE(s$paired)) " for paired samples" else "", ", which compares the observed change in PSI with the changes ",
         "between replicates of events with similar expression", if (isTRUE(s$gene_correction)) "; p-values were corrected for the number of events per gene" else "",
         ". ", s$sign_convention, ". An event was called differentially spliced when ", s$significance_rule, ": ",
         s$n_significant, " of ", s$n_events_tested, " tested events (", s$n_up, " with more and ", s$n_down,
         " with less inclusion in ", test, ") in ", s$n_genes_significant, " genes."), "",
  if (!is.null(cmp) && identical(cmp$status, "done")) paste0(
    "Events were matched with rMATS events from the same samples by their splice junctions (", cmp$rmats_count_type,
    " counts): ", cmp$n_matched, " events were tested by both tools, the change in PSI agreed with Spearman's rho = ",
    cmp$spearman_dpsi_matched, ", and ", cmp$n_sig_both, " events were significant in both (", cmp$rule, ")."), "",
  paste0("Analysis code: RNA-seq code bank ", CODE_BANK_VERSION, ", SUPPA module (", params$suppa_command, "); R ",
         R.version$major, ".", R.version$minor, "."), "",
  "References", "",
  "- Trincado JL, Entizne JC, Hysenaj G, Singh B, Skalic M, Elliott DJ, Eyras E (2018). SUPPA2: fast, accurate, and uncertainty-aware differential splicing analysis across multiple conditions. Genome Biology 19:40.",
  "- Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods 14:417-419.",
  if (!is.null(cmp) && identical(cmp$status, "done"))
    "- Shen S, Park JW, Lu ZX, Lin L, Henry MD, Wu YN, Zhou Q, Xing Y (2014). rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. PNAS 111:E5593-E5601."
)
writeLines(methods[!is.na(methods)], file.path(out_dir, "methods.md"))
write_step_provenance(run_dir, STEP, started, params_used = list(), outputs = list.files(out_dir))
log_info("Step SU05 done: methods.md written.")
