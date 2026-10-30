#!/usr/bin/env Rscript
# =============================================================================
# splicing/tests/run_splicing_tests.R  --  end-to-end tests of the splicing module
# =============================================================================
# Run from the repository root:  Rscript splicing/tests/run_splicing_tests.R
#
# Synthetic rMATS tables (all five event types, rMATS 4 column layout) with
# planted events, so the right answer is known:
#   A. a full run finishes and writes every expected file
#   B. exactly the planted events are significant, with the right direction
#   C. low-coverage events are not called significant even with a tiny FDR
#   D. enrichment finds the planted pathway (GMT gene sets, offline)
#   E. a rerun from saved parameters gives identical results; a stricter cutoff gives fewer events
#   F. bad inputs fail with the right check (no tables, --statoff output, missing columns)
#   G. with a DE run: genes are sorted into splicing only / expression only / both;
#      enrichment also tests genes that only gain or only lose inclusion
# Exit code 0 = all passed. Test runs go to splicing/tests/output/ (ignored by git).
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "splicing/tests/x"), "..", ".."))
})
setwd(REPO_ROOT)
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "splicing", "R", "splicing_utils.R"))

out_root <- file.path(REPO_ROOT, "splicing", "tests", "output")
unlink(out_root, recursive = TRUE)
dir.create(out_root, recursive = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")
results <- data.frame(test = character(0), status = character(0), detail = character(0))
record <- function(test, ok, detail = "") {
  results[nrow(results) + 1, ] <<- list(test, if (isTRUE(ok)) "PASS" else "FAIL", detail)
  cat(sprintf("  [%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", test, if (nzchar(detail)) paste0("- ", detail) else ""))
}
run_splicing <- function(args) {
  before <- list.dirs(out_root, recursive = FALSE)
  log <- tempfile(fileext = ".log")
  code <- system2(rscript, c("splicing/scripts/run_splicing.R", "--runs_dir", shQuote(out_root), args), stdout = log, stderr = log)
  after <- setdiff(list.dirs(out_root, recursive = FALSE), before)
  list(code = code, run_dir = if (length(after)) after[1] else NA, log = log)
}

# ---- Synthetic rMATS output -----------------------------------------------------------
set.seed(7)
n_genes <- 600
genes <- sprintf("GENE%03d", seq_len(n_genes))
planted_up <- genes[1:24]      # more inclusion in KO
planted_down <- genes[25:40]   # less inclusion in KO
low_cov <- genes[41:45]        # big change, tiny FDR, but too few reads
make_table <- function(et, n_events = 300) {
  g <- c(planted_up, planted_down, low_cov, sample(genes[46:n_genes], n_events - 45))
  base <- sample(1e5:9e7, n_events)
  coords <- switch(et,
    SE = data.frame(exonStart_0base = base + 500, exonEnd = base + 600, upstreamES = base, upstreamEE = base + 100,
                    downstreamES = base + 1000, downstreamEE = base + 1100),
    MXE = data.frame(`1stExonStart_0base` = base + 400, `1stExonEnd` = base + 500, `2ndExonStart_0base` = base + 700,
                     `2ndExonEnd` = base + 800, upstreamES = base, upstreamEE = base + 100, downstreamES = base + 1000,
                     downstreamEE = base + 1100, check.names = FALSE),
    A3SS = , A5SS = data.frame(longExonStart_0base = base, longExonEnd = base + 300, shortES = base + 100, shortEE = base + 300,
                               flankingES = base + 1000, flankingEE = base + 1100),
    RI = data.frame(riExonStart_0base = base, riExonEnd = base + 900, upstreamES = base, upstreamEE = base + 100,
                    downstreamES = base + 800, downstreamEE = base + 900))
  psi2 <- runif(n_events, 0.2, 0.8)
  dpsi <- rnorm(n_events, 0, 0.02)
  dpsi[g %in% planted_up] <- 0.3
  dpsi[g %in% planted_down] <- -0.3
  dpsi[g %in% low_cov] <- 0.4
  psi1 <- pmin(pmax(psi2 + dpsi, 0), 1)
  reads <- ifelse(g %in% low_cov, 3, 60)
  rep_vals <- function(p, r) vapply(seq_along(p), function(i) paste(round(pmin(pmax(p[i] + rnorm(3, 0, 0.02), 0), 1), 3), collapse = ","), "")
  ijc <- function(p, r) vapply(seq_along(p), function(i) paste(rpois(3, r[i] * p[i]), collapse = ","), "")
  sjc <- function(p, r) vapply(seq_along(p), function(i) paste(rpois(3, r[i] * (1 - p[i])), collapse = ","), "")
  pv <- runif(n_events)
  pv[g %in% c(planted_up, planted_down, low_cov)] <- 1e-9
  fdr <- p.adjust(pv, "BH")
  data.frame(ID = seq_len(n_events) - 1, GeneID = sprintf("\"ENSG%011d\"", match(g, genes)), geneSymbol = sprintf("\"%s\"", g),
             chr = "chr1", strand = sample(c("+", "-"), n_events, TRUE), coords, ID = seq_len(n_events) - 1,
             IJC_SAMPLE_1 = ijc(psi1, reads), SJC_SAMPLE_1 = sjc(psi1, reads), IJC_SAMPLE_2 = ijc(psi2, reads),
             SJC_SAMPLE_2 = sjc(psi2, reads), IncFormLen = 249, SkipFormLen = 99, PValue = signif(pv, 4), FDR = signif(fdr, 4),
             IncLevel1 = rep_vals(psi1), IncLevel2 = rep_vals(psi2), IncLevelDifference = round(dpsi, 3), check.names = FALSE)
}
rm_dir <- file.path(out_root, "rmats_output")
dir.create(rm_dir)
for (et in EVENT_TYPES) {
  tab <- make_table(et)
  # SE carries all planted events; other types get only background (no planted genes)
  if (et != "SE") {
    planted <- grepl("GENE0[0-3][0-9]\"|GENE04[0-5]\"", tab$geneSymbol)
    tab$IncLevelDifference[planted] <- 0; tab$PValue[planted] <- 0.9; tab$FDR[planted] <- 1
  }
  data.table::fwrite(tab, file.path(rm_dir, sprintf("%s.MATS.JCEC.txt", et)), sep = "\t", quote = FALSE)
}
gmt <- file.path(out_root, "test_sets.gmt")
writeLines(c(paste(c("PLANTED_SPLICING_SET", "test", genes[1:40], genes[200:215]), collapse = "\t"),
             paste(c("BACKGROUND_SET", "test", genes[300:360]), collapse = "\t")), gmt)
common <- c("--rmats_dir", shQuote(rm_dir), "--group_1_name", "KO", "--group_2_name", "WT",
            "--set", shQuote(sprintf("gene_set_source=gmt;gmt_files=[%s];group_1_samples=[KO1, KO2, KO3];group_2_samples=[WT1, WT2, WT3]", gmt)))

# ---- A. full run -------------------------------------------------------------------------
cat("\nA. full splicing run\n")
a <- run_splicing(c(common, "--run_name", "test_splicing"))
record("run exit code 0", a$code == 0, a$run_dir)
rd <- a$run_dir
expected <- c("S01_validate_rmats/validation_report.json", "S01_validate_rmats/events_validated.csv",
              "S02_splicing_results/splicing_events_all.csv", "S02_splicing_results/splicing_events_significant.csv",
              "S02_splicing_results/splicing_summary.json", "S02_splicing_results/event_counts.png",
              "S02_splicing_results/dpsi_volcano.png", "S02_splicing_results/psi_scatter.png",
              "S02_splicing_results/top_events_heatmap.png", "S03_splicing_enrichment/enrichment_summary.json",
              "S04_splicing_report/methods.md", "provenance.json", "run_summary.json", "params.yaml",
              "inputs/rmats/SE.MATS.JCEC.txt")
missing <- expected[!file.exists(file.path(rd, expected))]
record("every expected file exists", !length(missing), paste(missing, collapse = ", "))

# ---- B. planted events -----------------------------------------------------------------------
cat("\nB. planted events\n")
sig <- read.csv(file.path(rd, "S02_splicing_results", "splicing_events_significant.csv"))
record("all planted SE events found", setequal(sig$gene_symbol, c(planted_up, planted_down)),
       sprintf("%d significant", nrow(sig)))
record("directions right", all(sig$direction[sig$gene_symbol %in% planted_up] == "up") &&
         all(sig$direction[sig$gene_symbol %in% planted_down] == "down"))
record("no background events called", all(sig$event_type == "SE"))
ss <- read_json_file(file.path(rd, "S02_splicing_results", "splicing_summary.json"))
record("summary counts", ss$n_up == length(planted_up) && ss$n_down == length(planted_down),
       sprintf("up %d, down %d", ss$n_up, ss$n_down))

# ---- C. coverage filter --------------------------------------------------------------------------
cat("\nC. coverage filter\n")
all_ev <- read.csv(file.path(rd, "S02_splicing_results", "splicing_events_all.csv"))
lc <- all_ev[all_ev$gene_symbol %in% low_cov & all_ev$event_type == "SE", ]
record("low-coverage events not significant", nrow(lc) == length(low_cov) && !any(lc$significant) && all(!lc$passes_coverage))

# ---- D. enrichment ---------------------------------------------------------------------------------
cat("\nD. pathway enrichment\n")
ora <- read.csv(file.path(rd, "S03_splicing_enrichment", "test_sets", "ora_all_events.csv"))
record("planted pathway enriched", ora$pathway[1] == "PLANTED_SPLICING_SET" && ora$padj[1] < 0.05,
       sprintf("padj %.2g", ora$padj[1]))
record("background pathway not enriched", ora$padj[ora$pathway == "BACKGROUND_SET"] > 0.05)
record("methods text written", any(grepl("rMATS", readLines(file.path(rd, "S04_splicing_report", "methods.md")))))

# ---- E. reruns -----------------------------------------------------------------------------------------
cat("\nE. reruns\n")
e <- run_splicing(c("--based_on", shQuote(rd)))
e_all <- read.csv(file.path(e$run_dir, "S02_splicing_results", "splicing_events_all.csv"))
record("rerun gives identical events", e$code == 0 && isTRUE(all.equal(all_ev, e_all)))
s2 <- run_splicing(c("--based_on", shQuote(rd), "--set", shQuote("delta_psi_cutoff=0.35")))
s2s <- read_json_file(file.path(s2$run_dir, "S02_splicing_results", "splicing_summary.json"))
record("stricter dPSI cutoff gives fewer events", s2$code == 0 && s2s$n_significant < ss$n_significant,
       sprintf("%d < %d", s2s$n_significant, ss$n_significant))

# ---- F. bad inputs -------------------------------------------------------------------------------------
cat("\nF. bad inputs\n")
expect_fail <- function(label, check_prefix, dir, extra = character(0)) {
  r <- run_splicing(c("--rmats_dir", shQuote(dir), "--group_1_name", "KO", "--group_2_name", "WT", extra))
  rep <- read_json_file(file.path(r$run_dir, "S01_validate_rmats", "validation_report.json"))
  failed <- if (is.null(rep)) character(0) else rep$checks$id[rep$checks$status == "fail"]
  record(label, r$code == 1 && any(startsWith(failed, check_prefix)), paste(failed, collapse = ", "))
}
empty <- file.path(out_root, "empty"); dir.create(empty)
expect_fail("folder without rMATS tables", "files_found", empty)
statoff <- file.path(out_root, "statoff"); dir.create(statoff)
t <- read_rmats_table(file.path(rm_dir, "SE.MATS.JCEC.txt")); t$PValue <- NA; t$FDR <- NA
names(t)[names(t) == "ID.1"] <- "ID"
data.table::fwrite(t, file.path(statoff, "SE.MATS.JCEC.txt"), sep = "\t", quote = FALSE)
expect_fail("rMATS run with --statoff", "statistics_SE", statoff, c("--set", shQuote("event_types=[SE]")))
cols <- file.path(out_root, "cols"); dir.create(cols)
t2 <- read_rmats_table(file.path(rm_dir, "SE.MATS.JCEC.txt")); t2$IncLevelDifference <- NULL
data.table::fwrite(t2, file.path(cols, "SE.MATS.JCEC.txt"), sep = "\t", quote = FALSE)
expect_fail("missing columns", "columns_SE", cols, c("--set", shQuote("event_types=[SE]")))
expect_fail("JC tables asked for, only JCEC present", "files_found", rm_dir, c("--set", shQuote("count_type=JC")))

# ---- G. splicing vs expression, directional gene lists ---------------------------------------------------
cat("\nG. splicing vs expression\n")
de_dir <- file.path(out_root, "fake_de_run", "04_differential_expression"); dir.create(de_dir, recursive = TRUE)
de_genes <- c(genes[1:12], genes[300:320])          # half the 'up' splicing genes also change expression; plus DE-only genes
data.table::fwrite(data.frame(gene_id = sprintf("ENSG%011d", seq_len(n_genes)), gene_symbol = genes,
                              log2FoldChange = ifelse(genes %in% de_genes, 2, 0.05), padj = ifelse(genes %in% de_genes, 1e-6, 0.8),
                              significant = genes %in% de_genes), file.path(de_dir, "de_results_all.csv"))
g <- run_splicing(c(common, "--set", shQuote(sprintf("de_run=%s;gene_set_source=gmt;gmt_files=[%s]", dirname(de_dir), gmt))))
sve <- read.csv(file.path(g$run_dir, "S02_splicing_results", "splicing_vs_expression.csv"))
cat_of <- function(x) sve$category[sve$gene_symbol == x]
record("splicing vs expression written", g$code == 0 && file.exists(file.path(g$run_dir, "S02_splicing_results", "splicing_vs_expression.png")))
record("categories right", cat_of("GENE001") == "splicing and expression" && cat_of("GENE020") == "splicing only" &&
         cat_of("GENE300") == "expression only", paste(cat_of("GENE001"), cat_of("GENE020"), cat_of("GENE300"), sep = " / "))
record("dPSI distribution figure", file.exists(file.path(g$run_dir, "S02_splicing_results", "dpsi_distribution.png")))
en <- read_json_file(file.path(g$run_dir, "S03_splicing_enrichment", "enrichment_summary.json"))
record("directional gene lists tested", en$n_genes$more_inclusion == length(planted_up) && en$n_genes$less_inclusion == length(planted_down) &&
         file.exists(file.path(g$run_dir, "S03_splicing_enrichment", "test_sets", "ora_more_inclusion.csv")),
       sprintf("%s more, %s less", en$n_genes$more_inclusion, en$n_genes$less_inclusion))

cat("\n================================================================\n")
n_fail <- sum(results$status == "FAIL")
cat(sprintf(" %d tests: %d passed, %d failed\n", nrow(results), sum(results$status == "PASS"), n_fail))
if (n_fail) print(results[results$status == "FAIL", ], row.names = FALSE)
cat("================================================================\n")
write.csv(results, file.path(out_root, "test_results.csv"), row.names = FALSE)
quit(save = "no", status = if (n_fail) 1 else 0)
