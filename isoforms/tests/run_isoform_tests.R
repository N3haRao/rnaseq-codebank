#!/usr/bin/env Rscript
# =============================================================================
# isoforms/tests/run_isoform_tests.R  --  end-to-end tests of the isoform module
# =============================================================================
# Run from the repository root:  Rscript isoforms/tests/run_isoform_tests.R
# Needs DESeq2, DEXSeq, stageR and tximport. Output: isoforms/tests/output/.
#
# Synthetic Salmon output (quant.sf per sample) for 800 genes with 1-4
# transcripts each, 4 WT and 4 KO samples, with two kinds of planted change:
#   expression genes: every transcript of the gene 4x higher in KO (shares unchanged)
#   switch genes:     the main transcript goes from 80% to 20% of the gene in KO,
#                     the gene total unchanged
# Expected: DEI finds the expression genes' transcripts (and the switching
# transcripts); DUI finds the switch genes and NOT the expression-only genes.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "isoforms/tests/x"), "..", ".."))
})
setwd(REPO_ROOT)
source(file.path(REPO_ROOT, "R", "utils.R"))
missing_pkgs <- c("DESeq2", "DEXSeq", "stageR", "tximport")[!vapply(c("DESeq2", "DEXSeq", "stageR", "tximport"), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) { cat("SKIP: needs", paste(missing_pkgs, collapse = ", "), "\n"); quit(save = "no", status = 0) }

out <- file.path(REPO_ROOT, "isoforms", "tests", "output")
unlink(out, recursive = TRUE)
dir.create(file.path(out, "salmon"), recursive = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")
results <- data.frame(test = character(0), status = character(0), detail = character(0))
record <- function(test, ok, detail = "") {
  results[nrow(results) + 1, ] <<- list(test, if (isTRUE(ok)) "PASS" else "FAIL", detail)
  cat(sprintf("  [%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", test, if (nzchar(detail)) paste0("- ", detail) else ""))
}
run_iso <- function(args) {
  before <- list.dirs(out, recursive = FALSE)
  log <- tempfile(fileext = ".log")
  code <- system2(rscript, c("isoforms/scripts/run_isoforms.R", "--runs_dir", shQuote(out), args), stdout = log, stderr = log)
  if (code != 0) cat(paste0("    | ", tail(readLines(log), 8), "\n"), sep = "")
  list(code = code, run_dir = setdiff(list.dirs(out, recursive = FALSE), before)[1])
}

# ---- Synthetic transcriptome and Salmon output -------------------------------------------------
set.seed(3)
n_genes <- 800
n_tx <- sample(1:4, n_genes, TRUE, prob = c(0.3, 0.35, 0.2, 0.15))
genes <- sprintf("ENSG%011d", seq_len(n_genes))
tx <- data.frame(gene_id = rep(genes, n_tx), gene_name = rep(sprintf("GENE%03d", seq_len(n_genes)), n_tx))
tx$transcript_id <- sprintf("ENST%011d.1", seq_len(nrow(tx)))
tx$length <- sample(800:4000, nrow(tx), TRUE)
multi <- genes[n_tx >= 2]
expr_genes <- multi[1:25]            # all transcripts 4x up in KO
switch_genes <- multi[26:50]         # main transcript 0.8 -> 0.2 share
gene_mean <- setNames(round(exp(rnorm(n_genes, log(800), 0.8))), genes)
share_wt <- unlist(lapply(n_tx, function(k) { if (k == 1) 1 else { p <- c(0.8, rep(0.2 / (k - 1), k - 1)); p } }))
tx$share_wt <- share_wt
tx$share_ko <- tx$share_wt
for (g in switch_genes) {
  i <- which(tx$gene_id == g); k <- length(i)
  tx$share_ko[i] <- c(0.2, rep(0.8 / (k - 1), k - 1))
}
samples <- data.frame(sample_id = c(paste0("WT", 1:4), paste0("KO", 1:4)), genotype = rep(c("WT", "KO"), each = 4),
                      batch = rep(c("b1", "b2"), 4))
for (i in seq_len(nrow(samples))) {
  ko <- samples$genotype[i] == "KO"
  mu <- gene_mean[tx$gene_id] * (if (ko) tx$share_ko else tx$share_wt) * ifelse(ko & tx$gene_id %in% expr_genes, 4, 1)
  reads <- rnbinom(nrow(tx), mu = mu, size = 20)
  eff <- tx$length - 200
  tpm <- reads / eff; tpm <- 1e6 * tpm / sum(tpm)
  d <- file.path(out, "salmon", samples$sample_id[i]); dir.create(d)
  write.table(data.frame(Name = tx$transcript_id, Length = tx$length, EffectiveLength = eff, TPM = round(tpm, 4), NumReads = reads),
              file.path(d, "quant.sf"), sep = "\t", quote = FALSE, row.names = FALSE)
}
write.csv(samples, file.path(out, "samples.csv"), row.names = FALSE)
gtf <- c(sprintf('chr1\ttest\ttranscript\t%d\t%d\t.\t+\t.\tgene_id "%s"; transcript_id "%s"; gene_name "%s";',
                 seq_len(nrow(tx)) * 10, seq_len(nrow(tx)) * 10 + 5, tx$gene_id, tx$transcript_id, tx$gene_name))
writeLines(gtf, file.path(out, "genes.gtf"))
common <- c("--quant_dir", shQuote(file.path(out, "salmon")), "--metadata_file", shQuote(file.path(out, "samples.csv")),
            "--gtf_file", shQuote(file.path(out, "genes.gtf")), "--condition", "genotype", "--reference_level", "WT", "--test_level", "KO")

# ---- A. full run -------------------------------------------------------------------------------------
cat("\nA. full isoform run\n")
a <- run_iso(c(common, "--run_name", "test_isoforms", "--set", shQuote("covariates=[batch]")))
rd <- a$run_dir
record("run exit code 0", a$code == 0, rd)
expected <- c("I01_validate_isoforms/validation_report.json", "I01_validate_isoforms/txi_counts.rds",
              "I02_differential_isoform_expression/dei_results_all.csv", "I02_differential_isoform_expression/dei_volcano.png",
              "I03_differential_isoform_usage/dui_transcripts_all.csv", "I03_differential_isoform_usage/dui_genes.csv",
              "I03_differential_isoform_usage/dui_volcano.png", "I03_differential_isoform_usage/dui_top_genes_usage.png",
              "I04_isoform_report/methods.md", "run_summary.json", "provenance.json", "inputs/quant/KO1/quant.sf")
miss <- expected[!file.exists(file.path(rd, expected))]
record("every expected file exists", !length(miss), paste(miss, collapse = ", "))

# ---- B. DEI ----------------------------------------------------------------------------------------------
cat("\nB. differential isoform expression\n")
dei <- read.csv(file.path(rd, "I02_differential_isoform_expression", "dei_results_all.csv"))
de_genes <- unique(dei$gene_id[dei$significant])
expr_hit <- mean(expr_genes %in% de_genes)
record("expression genes found by DEI", expr_hit >= 0.9, sprintf("%.0f%% of planted", 100 * expr_hit))
up_ok <- all(dei$direction[dei$significant & dei$gene_id %in% expr_genes] == "up")
record("their transcripts are higher in KO", up_ok)
false_dei <- setdiff(de_genes, c(expr_genes, switch_genes))
record("few false DEI genes", length(false_dei) <= 3, sprintf("%d unplanted genes called", length(false_dei)))

# ---- C. DUI -------------------------------------------------------------------------------------------------
cat("\nC. differential isoform usage\n")
dg <- read.csv(file.path(rd, "I03_differential_isoform_usage", "dui_genes.csv"))
dui_genes <- dg$gene_id[dg$significant]
sw_hit <- mean(switch_genes %in% dui_genes)
record("switch genes found by DUI", sw_hit >= 0.9, sprintf("%.0f%% of planted", 100 * sw_hit))
record("expression-only genes NOT called by DUI", !any(expr_genes %in% dui_genes),
       paste(intersect(expr_genes, dui_genes), collapse = ", "))
record("few false DUI genes", length(setdiff(dui_genes, switch_genes)) <= 2, sprintf("%d unplanted genes", length(setdiff(dui_genes, switch_genes))))
dt <- read.csv(file.path(rd, "I03_differential_isoform_usage", "dui_transcripts_all.csv"))
main_tx <- tx$transcript_id[match(switch_genes, tx$gene_id)]
mt <- dt[dt$transcript_id %in% main_tx & dt$significant, ]
record("main transcript of switch genes loses share (~ -0.6)", nrow(mt) >= 20 && all(mt$direction == "down") && abs(median(mt$delta_usage) + 0.6) < 0.1,
       sprintf("%d transcripts, median delta %.2f", nrow(mt), median(mt$delta_usage)))
record("isoform switches reported", sum(dg$isoform_switch) >= 20, sprintf("%d switches", sum(dg$isoform_switch)))
record("methods text mentions DEXSeq and stageR", all(c("DEXSeq", "stageR") %in% unlist(strsplit(paste(readLines(file.path(rd, "I04_isoform_report", "methods.md")), collapse = " "), "[^A-Za-z]+"))))

# ---- D. rerun and bad inputs ------------------------------------------------------------------------------------
cat("\nD. rerun and bad inputs\n")
d <- run_iso(c("--based_on", shQuote(rd), "--steps", "1-2"))
dei2 <- read.csv(file.path(d$run_dir, "I02_differential_isoform_expression", "dei_results_all.csv"))
record("rerun gives identical DEI", d$code == 0 && isTRUE(all.equal(dei, dei2)))
bad <- samples; bad$sample_id[1] <- "WTX"
write.csv(bad, file.path(out, "bad_samples.csv"), row.names = FALSE)
b1 <- run_iso(c(sub("samples.csv", "bad_samples.csv", common, fixed = TRUE)))
rep1 <- read_json_file(file.path(b1$run_dir, "I01_validate_isoforms", "validation_report.json"))
record("missing quant.sf for a sample fails", b1$code == 1 && "quant_files" %in% rep1$checks$id[rep1$checks$status == "fail"])
conf <- samples; conf$batch <- conf$genotype
write.csv(conf, file.path(out, "conf_samples.csv"), row.names = FALSE)
b2 <- run_iso(c(sub("samples.csv", "conf_samples.csv", common, fixed = TRUE), "--set", shQuote("covariates=[batch]")))
rep2 <- read_json_file(file.path(b2$run_dir, "I01_validate_isoforms", "validation_report.json"))
record("covariate confounded with the groups fails", b2$code == 1 && "design_full_rank" %in% rep2$checks$id[rep2$checks$status == "fail"])
writeLines(sub("ENST", "XXXX", gtf), file.path(out, "wrong.gtf"))
b3 <- run_iso(c(sub("genes.gtf", "wrong.gtf", common, fixed = TRUE)))
rep3 <- read_json_file(file.path(b3$run_dir, "I01_validate_isoforms", "validation_report.json"))
record("annotation that does not match fails", b3$code == 1 && "tx2gene" %in% rep3$checks$id[rep3$checks$status == "fail"])

cat("\n================================================================\n")
n_fail <- sum(results$status == "FAIL")
cat(sprintf(" %d tests: %d passed, %d failed\n", nrow(results), sum(results$status == "PASS"), n_fail))
if (n_fail) print(results[results$status == "FAIL", ], row.names = FALSE)
cat("================================================================\n")
quit(save = "no", status = if (n_fail) 1 else 0)
