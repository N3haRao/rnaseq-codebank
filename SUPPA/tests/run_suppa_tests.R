#!/usr/bin/env Rscript
# =============================================================================
# suppa/tests/run_suppa_tests.R  --  end-to-end tests of the SUPPA2 module
# =============================================================================
# Run from the repository root:  Rscript suppa/tests/run_suppa_tests.R
# Needs SUPPA2 (suppa/install/install_suppa.sh, or SUPPA_COMMAND=/path/to/suppa.py);
# skipped without it. Output: suppa/tests/output/.
#
# Synthetic annotation and Salmon output, 3 WT and 3 KO samples:
#   300 genes with a cassette exon (skipped exon, SE): the first 30 gain inclusion in KO (+0.4)
#   100 genes with a retained intron (RI): the first 10 lose retention in KO (-0.35)
# Expected: SUPPA2 finds the planted events with the right sign, few others; the rMATS
# comparison matches the SE events by their junctions; bad inputs fail with clear checks.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "suppa/tests/x"), "..", ".."))
})
setwd(REPO_ROOT)
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "suppa", "R", "suppa_utils.R"))
cmd <- find_suppa(list())
if (is.null(cmd)) { cat("SKIP: SUPPA2 not found (suppa/install/install_suppa.sh, or set SUPPA_COMMAND)\n"); quit(save = "no", status = 0) }
cat("SUPPA2:", paste(cmd, collapse = " "), "\n")

out <- file.path(REPO_ROOT, "suppa", "tests", "output")
unlink(out, recursive = TRUE)
dir.create(file.path(out, "salmon"), recursive = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")
results <- data.frame(test = character(0), status = character(0), detail = character(0))
record <- function(test, ok, detail = "") {
  if (!length(detail) || is.na(detail[1])) detail <- "" else detail <- as.character(detail[1])
  results[nrow(results) + 1, ] <<- list(test, if (isTRUE(ok)) "PASS" else "FAIL", detail)
  cat(sprintf("  [%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", test, if (nzchar(detail)) paste0("- ", detail) else ""))
}
run_suppa_module <- function(args, env = character(0)) {
  before <- list.dirs(out, recursive = FALSE)
  log <- tempfile(fileext = ".log")
  code <- system2(rscript, c("suppa/scripts/run_suppa.R", "--runs_dir", shQuote(out), args), stdout = log, stderr = log, env = env)
  if (code != 0) cat(paste0("    | ", tail(readLines(log), 8), "\n"), sep = "")
  list(code = code, run_dir = setdiff(list.dirs(out, recursive = FALSE), before)[1], log = log)
}

# ---- Synthetic annotation ---------------------------------------------------------------------
set.seed(11)
gtf <- character(0); se_genes <- character(0); ri_genes <- character(0); tx_rows <- list()
add_tx <- function(gid, gname, tid, exons, chr, start) {
  s <- start + exons[, 1]; e <- start + exons[, 2]
  gtf <<- c(gtf, sprintf('%s\tT\ttranscript\t%d\t%d\t.\t+\t.\tgene_id "%s"; transcript_id "%s"; gene_name "%s";', chr, min(s), max(e), gid, tid, gname),
            sprintf('%s\tT\texon\t%d\t%d\t.\t+\t.\tgene_id "%s"; transcript_id "%s"; gene_name "%s";', chr, s, e, gid, tid, gname))
}
k <- 0
for (i in 1:400) {
  gid <- sprintf("ENSG%011d.1", i); gname <- sprintf("GENE%03d", i); start <- i * 10000
  gtf <- c(gtf, sprintf('chr1\tT\tgene\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_name "%s";', start + 100, start + 1000, gid, gname))
  t1 <- sprintf("ENST%011d.1", 2 * i - 1); t2 <- sprintf("ENST%011d.1", 2 * i)
  if (i <= 300) {          # cassette exon: inclusion form (3 exons) and skipping form (2 exons)
    add_tx(gid, gname, t1, rbind(c(100, 200), c(400, 500), c(800, 1000)), "chr1", start)
    add_tx(gid, gname, t2, rbind(c(100, 200), c(800, 1000)), "chr1", start)
    se_genes <- c(se_genes, gid)
  } else {                 # retained intron: one exon 100-1000 vs exons 100-200 + 800-1000
    add_tx(gid, gname, t1, rbind(c(100, 1000)), "chr1", start)
    add_tx(gid, gname, t2, rbind(c(100, 200), c(800, 1000)), "chr1", start)
    ri_genes <- c(ri_genes, gid)
  }
  tx_rows[[i]] <- data.frame(gene_id = gid, gene_name = gname, inc = t1, alt = t2, start = start)
}
writeLines(gtf, file.path(out, "genes.gtf"))
txg <- do.call(rbind, tx_rows)

# ---- Synthetic Salmon output: PSI = share of the inclusion (t1) form --------------------------------
samples <- c(paste0("WT_", 1:3), paste0("KO_", 1:3))
base_psi <- runif(400, 0.25, 0.65)
shift <- rep(0, 400); shift[1:30] <- 0.4; shift[301:310] <- -0.35
gene_tpm <- exp(rnorm(400, log(60), 0.5))
psi_mat <- matrix(NA_real_, 400, length(samples), dimnames = list(NULL, samples))
for (s in samples) {
  ko <- startsWith(s, "KO")
  psi <- pmin(0.98, pmax(0.02, base_psi + (if (ko) shift else 0) + rnorm(400, 0, 0.03)))
  psi_mat[, s] <- psi
  tot <- gene_tpm * exp(rnorm(400, 0, 0.1))
  tpm <- c(rbind(tot * psi, tot * (1 - psi)))
  ids <- c(rbind(txg$inc, txg$alt)); genes <- rep(txg$gene_id, each = 2)
  dir.create(file.path(out, "salmon", s))
  data.table::fwrite(data.frame(Name = paste0(ids, "|", genes, "|-|-|x-201|x|1000|protein_coding|"), Length = 1000,
                                EffectiveLength = 850, TPM = round(tpm, 4), NumReads = round(tpm * 20)),
                     file.path(out, "salmon", s, "quant.sf"), sep = "\t")
}
write.csv(data.frame(sample = samples, genotype = rep(c("WT", "KO"), each = 3)), file.path(out, "samples.csv"), row.names = FALSE)

# ---- Synthetic rMATS SE table for the same comparison (--b1 = KO) ----------------------------------
se <- txg[1:300, ]
fdr <- ifelse(seq_len(300) <= 30, 1e-6, runif(300, 0.2, 1))
data.table::fwrite(data.frame(ID = 1:300, GeneID = se$gene_id, geneSymbol = se$gene_name, chr = "chr1", strand = "+",
                              exonStart_0base = se$start + 399, exonEnd = se$start + 500, upstreamES = se$start + 99,
                              upstreamEE = se$start + 200, downstreamES = se$start + 799, downstreamEE = se$start + 1000,
                              PValue = fdr, FDR = fdr, IncLevel1 = "0.5,0.5,0.5", IncLevel2 = "0.5,0.5,0.5",
                              IncLevelDifference = round(rowMeans(psi_mat[1:300, 4:6]) - rowMeans(psi_mat[1:300, 1:3]) + rnorm(300, 0, 0.01), 3)),
                   file.path(out, "SE.MATS.JCEC.txt"), sep = "\t")
dir.create(file.path(out, "rmats")); invisible(file.rename(file.path(out, "SE.MATS.JCEC.txt"), file.path(out, "rmats", "SE.MATS.JCEC.txt")))

# ---- 1. Full run ---------------------------------------------------------------------------------
common <- c("--quant_dir", shQuote(file.path(out, "salmon")), "--metadata_file", shQuote(file.path(out, "samples.csv")),
            "--gtf_file", shQuote(file.path(out, "genes.gtf")), "--condition", "genotype", "--reference_level", "WT", "--test_level", "KO")
cat("\nFull run (SE + RI, rMATS comparison):\n")
r <- run_suppa_module(c(common, "--rmats_dir", shQuote(file.path(out, "rmats")), "--set", shQuote("event_types=[SE,RI]")))
summ <- read_json_file(file.path(r$run_dir, "run_summary.json"))
record("run finishes with status success", r$code == 0 && identical(summ$status, "success"), summ$status %||% "")
evs <- read_json_file(file.path(r$run_dir, "SU02_suppa_events_psi", "events_summary.json"))
record("events found in the annotation (300 SE, 100 RI)", identical(as.integer(evs$events_by_type$SE), 300L) && identical(as.integer(evs$events_by_type$RI), 100L),
       sprintf("SE %s, RI %s", evs$events_by_type$SE, evs$events_by_type$RI))
all <- data.table::fread(file.path(r$run_dir, "SU03_suppa_diffsplice", "suppa_events_all.csv"), data.table = FALSE)
sig <- all[all$significant, ]
planted_se <- txg$gene_id[1:30]; planted_ri <- txg$gene_id[301:310]
hit_se <- sig$gene_id %in% planted_se & sig$event_type == "SE"
hit_ri <- sig$gene_id %in% planted_ri & sig$event_type == "RI"
record("planted SE events found (>= 80%)", sum(hit_se) >= 24, sprintf("%d of 30", sum(hit_se)))
record("planted RI events found (>= 70%)", sum(hit_ri) >= 7, sprintf("%d of 10", sum(hit_ri)))
record("sign: more inclusion in KO = positive dPSI", all(sig$dpsi[hit_se] > 0) && all(sig$dpsi[hit_ri] < 0))
fp <- sum(!(hit_se | hit_ri))
record("few false positives (<= 5% of the 360 unchanged events)", fp <= 18, sprintf("%d", fp))
record("gene names come from the GTF", all(grepl("^GENE", sig$gene_name)))
psum <- read_json_file(file.path(r$run_dir, "SU03_suppa_diffsplice", "suppa_summary.json"))
record("summary counts agree with the table", identical(as.integer(psum$n_significant), nrow(sig)) && psum$n_up + psum$n_down == nrow(sig))
for (f in c("event_counts.png", "dpsi_volcano.png", "psi_scatter.png"))
  record(paste("figure", f), file.exists(file.path(r$run_dir, "SU03_suppa_diffsplice", f)))
cmp <- read_json_file(file.path(r$run_dir, "SU04_suppa_vs_rmats", "comparison_summary.json"))
record("rMATS comparison matches the SE events by junctions", identical(cmp$status, "done") && cmp$n_matched >= 290, sprintf("%s matched", cmp$n_matched))
record("rMATS and SUPPA2 agree (Spearman > 0.7; >= 24 significant in both)", isTRUE(cmp$spearman_dpsi_matched > 0.7) && cmp$n_sig_both >= 24,
       sprintf("rho %s, %s in both", cmp$spearman_dpsi_matched, cmp$n_sig_both))
record("methods text written", file.exists(file.path(r$run_dir, "SU05_suppa_report", "methods.md")) &&
         any(grepl("SUPPA2", readLines(file.path(r$run_dir, "SU05_suppa_report", "methods.md")))))

# ---- 2. Rerun from the run's parameters, with a stricter cutoff -------------------------------------
r2 <- run_suppa_module(c("--based_on", shQuote(r$run_dir), "--set", shQuote("dpsi_cutoff=0.3")))
s2 <- read_json_file(file.path(r2$run_dir, "SU03_suppa_diffsplice", "suppa_summary.json"))
record("--based_on rerun with dpsi_cutoff 0.3 finds fewer or equal events", r2$code == 0 && s2$n_significant <= psum$n_significant,
       sprintf("%s vs %s", s2$n_significant, psum$n_significant))

# ---- 3. Bad inputs ----------------------------------------------------------------------------------
cat("\nBad inputs:\n")
check_fail <- function(run, id) {
  rep <- read_json_file(file.path(run$run_dir, "SU01_validate_suppa", "validation_report.json"))
  run$code == 1 && id %in% rep$checks$id[rep$checks$status == "fail"]
}
dir.create(file.path(out, "salmon_missing")); for (s in samples[-6]) file.copy(file.path(out, "salmon", s), file.path(out, "salmon_missing"), recursive = TRUE)
b1 <- run_suppa_module(sub(file.path(out, "salmon"), file.path(out, "salmon_missing"), common, fixed = TRUE))
record("a sample without quant.sf fails", check_fail(b1, "quant_files"))
writeLines(gsub("ENST", "XXXX", gtf), file.path(out, "wrong.gtf"))
b2 <- run_suppa_module(sub("genes.gtf", "wrong.gtf", common, fixed = TRUE))
record("an annotation that does not match fails", check_fail(b2, "transcript_ids"))
b3 <- run_suppa_module(c(common, "--set", shQuote("suppa_command=/no/such/suppa.py")),
                       env = c("SUPPA_COMMAND=", "SUPPA_ENV=", "COPILOT_HOME=", paste0("PATH=", "/usr/bin:/bin")))
record("missing SUPPA2 fails with a clear message", check_fail(b3, "suppa"))
b4 <- run_suppa_module(replace(common, which(common == "WT"), "XX"))
record("an unknown reference group fails", check_fail(b4, "reference_level"))

cat("\n================================================================\n")
n_fail <- sum(results$status == "FAIL")
cat(sprintf(" %d tests: %d passed, %d failed\n", nrow(results), sum(results$status == "PASS"), n_fail))
if (n_fail) print(results[results$status == "FAIL", ], row.names = FALSE)
cat("================================================================\n")
quit(save = "no", status = if (n_fail) 1 else 0)
