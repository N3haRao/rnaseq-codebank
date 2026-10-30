#!/usr/bin/env Rscript
# =============================================================================
# later/tests/run_later_tests.R  --  end-to-end test of the FASTQ / BAM scripts
# =============================================================================
# Run from the repository root:  Rscript later/tests/run_later_tests.R
# Needs fastp, salmon, STAR, samtools and rmats.py on the PATH (e.g. one conda env:
#   conda create -n rnaseq-later -c conda-forge -c bioconda fastp salmon star samtools rmats)
# and the R package tximport. Missing tools: the test says so and stops (exit 0).
#
# A tiny synthetic genome with two genes. GENEA has a cassette exon (E2):
# WT samples include it in 80% of transcripts, KO samples in 20%. Paired-end,
# reverse-stranded (dUTP-like) reads are simulated, then:
#   U00 fastp keeps nearly all reads  U01 read length -> 75
#   U02 Salmon detects a reverse-stranded library
#   U03 gene counts for 2 genes      U06 STAR aligns; its counts also say "reverse"
#   U04 rMATS finds the E2 skipping
#   splicing module: the E2 event is significant with less inclusion in KO
#   U05 draws its sashimi plot
#   run_upstream.R does all of the above from one config file
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "later/tests/x"), "..", ".."))
})
setwd(REPO_ROOT)
source(file.path(REPO_ROOT, "R", "utils.R"))

tools <- c("fastp", "salmon", "STAR", "samtools", "rmats.py")
missing_tools <- tools[!nzchar(Sys.which(tools))]
if (length(missing_tools) || !requireNamespace("tximport", quietly = TRUE)) {
  cat("SKIP: needs", paste(c(missing_tools, if (!requireNamespace("tximport", quietly = TRUE)) "R package tximport"), collapse = ", "), "\n")
  quit(save = "no", status = 0)
}

out <- file.path(REPO_ROOT, "later", "tests", "output")
unlink(out, recursive = TRUE)
dir.create(file.path(out, "fastq"), recursive = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")
results <- data.frame(test = character(0), status = character(0), detail = character(0))
record <- function(test, ok, detail = "") {
  results[nrow(results) + 1, ] <<- list(test, if (isTRUE(ok)) "PASS" else "FAIL", detail)
  cat(sprintf("  [%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", test, if (nzchar(detail)) paste0("- ", detail) else ""))
}
run_r <- function(script, args) {
  log <- tempfile(fileext = ".log")
  code <- system2(rscript, c(script, args), stdout = log, stderr = log)
  if (code != 0) cat(paste0("    | ", tail(readLines(log), 12), "\n"), sep = "")
  code
}
sh <- function(cmd, args, log) system2(cmd, args, stdout = log, stderr = log)

# ---- Synthetic genome, annotation and reads ---------------------------------------------------
set.seed(11)
bases <- c("A", "C", "G", "T")
genome <- sample(bases, 12000, TRUE)
put <- function(pos, s) genome[pos:(pos + nchar(s) - 1)] <<- strsplit(s, "")[[1]]
exons <- list(E1 = c(1001, 1300), E2 = c(3001, 3150), E3 = c(5001, 5400), F1 = c(8001, 8600), F2 = c(9501, 10100))
# canonical GT...AG introns (on the + strand for GENEA; CT...AC on the genome for the - strand GENEB)
for (iv in list(c(1301, 3000), c(3151, 5000), c(1301, 5000))) { put(iv[1], "GT"); put(iv[2] - 1, "AG") }
put(8601, "CT"); put(9499, "AC")
revcomp <- function(s) paste(rev(c(A = "T", C = "G", G = "C", T = "A")[strsplit(s, "")[[1]]]), collapse = "")
seq_of <- function(ex) paste(genome[ex[1]:ex[2]], collapse = "")
tx_seq <- list(T1 = paste0(seq_of(exons$E1), seq_of(exons$E2), seq_of(exons$E3)),
               T2 = paste0(seq_of(exons$E1), seq_of(exons$E3)),
               T3 = revcomp(paste0(seq_of(exons$F1), seq_of(exons$F2))))
writeLines(c(">chr1", substring(paste(genome, collapse = ""), seq(1, 12000, 60), seq(60, 12000, 60))), file.path(out, "genome.fa"))
writeLines(unlist(lapply(names(tx_seq), function(t) c(paste0(">", t), tx_seq[[t]]))), file.path(out, "transcripts.fa"))
gtf_line <- function(feat, ex, strand, gid, gname, tid) sprintf('chr1\ttest\t%s\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; gene_name "%s";',
                                                            feat, ex[1], ex[2], strand, gid, tid, gname)
gtf <- c(gtf_line("gene", c(1001, 5400), "+", "ENSG00000000001", "GENEA", "")[0],
         gtf_line("transcript", c(1001, 5400), "+", "ENSG00000000001", "GENEA", "T1"),
         vapply(c("E1", "E2", "E3"), function(e) gtf_line("exon", exons[[e]], "+", "ENSG00000000001", "GENEA", "T1"), ""),
         gtf_line("transcript", c(1001, 5400), "+", "ENSG00000000001", "GENEA", "T2"),
         vapply(c("E1", "E3"), function(e) gtf_line("exon", exons[[e]], "+", "ENSG00000000001", "GENEA", "T2"), ""),
         gtf_line("transcript", c(8001, 10100), "-", "ENSG00000000002", "GENEB", "T3"),
         vapply(c("F2", "F1"), function(e) gtf_line("exon", exons[[e]], "-", "ENSG00000000002", "GENEB", "T3"), ""))
writeLines(gtf, file.path(out, "genes.gtf"))

simulate <- function(sample, p_t1, n = 4000, rl = 75) {
  tx <- sample(c("T1", "T2", "T3"), n, TRUE, prob = c(0.6 * p_t1, 0.6 * (1 - p_t1), 0.4))
  r1 <- r2 <- character(n)
  for (i in seq_len(n)) {
    s <- tx_seq[[tx[i]]]
    fl <- min(nchar(s), max(rl + 10, round(rnorm(1, 250, 25))))
    st <- sample.int(nchar(s) - fl + 1, 1)
    frag <- substr(s, st, st + fl - 1)
    r1[i] <- substr(revcomp(frag), 1, rl)    # reverse-stranded library: read 1 from the opposite strand
    r2[i] <- substr(frag, 1, rl)
  }
  q <- strrep("I", rl)
  for (k in 1:2) {
    con <- gzfile(file.path(out, "fastq", sprintf("%s_R%d.fastq.gz", sample, k)), "w")
    writeLines(as.vector(rbind(sprintf("@%s_%d/%d", sample, seq_len(n), k), if (k == 1) r1 else r2, "+", q)), con)
    close(con)
  }
}
samples <- data.frame(sample = c(paste0("WT", 1:3), paste0("KO", 1:3)), genotype = rep(c("WT", "KO"), each = 3))
for (i in seq_len(nrow(samples))) simulate(samples$sample[i], if (samples$genotype[i] == "WT") 0.8 else 0.2)
sheet <- do.call(rbind, lapply(seq_len(nrow(samples)), function(i) data.frame(
  file_name = sprintf("%s_R%d.fastq.gz", samples$sample[i], 1:2), sample = samples$sample[i],
  genotype = samples$genotype[i], read_orientation = c("R1", "R2"))))
write.table(sheet, file.path(out, "samples.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("Synthetic data written.\n")

# ---- U00 fastp ---------------------------------------------------------------------------------
cat("\nU00 QC and trimming\n")
code <- run_r("later/scripts/U00_fastq_qc_trim.R", c("--sample_sheet", file.path(out, "samples.tsv"), "--fastq_dir", file.path(out, "fastq"),
                                                     "--out_dir", file.path(out, "trimmed"), "--threads", "2"))
qc <- if (file.exists(file.path(out, "trimmed", "qc_summary.csv"))) read.csv(file.path(out, "trimmed", "qc_summary.csv")) else NULL
record("fastp ran for every sample", code == 0 && !is.null(qc) && nrow(qc) == 6 && all(qc$status == "success"))
record("nearly all reads kept", !is.null(qc) && all(qc$percent_kept > 95), if (!is.null(qc)) sprintf("min %.1f%%", min(qc$percent_kept)) else "")
ts <- if (file.exists(file.path(out, "trimmed", "samples_trimmed.tsv"))) read.delim(file.path(out, "trimmed", "samples_trimmed.tsv")) else NULL
record("trimmed sample sheet written", !is.null(ts) && nrow(ts) == 12 && all(file.exists(file.path(out, "trimmed", ts$file_name))) && "genotype" %in% names(ts))

# ---- U01 read length ----------------------------------------------------------------------------
cat("\nU01 read length\n")
code <- run_r("later/scripts/U01_read_length.R", c("--sample_sheet", file.path(out, "samples.tsv"), "--fastq_dir", file.path(out, "fastq"),
                                                   "--out_dir", file.path(out, "read_length")))
rl <- suppressWarnings(as.integer(readLines(file.path(out, "read_length", "read_length.txt"))))
record("read length 75", code == 0 && identical(rl, 75L), as.character(rl))

# ---- U02 Salmon ---------------------------------------------------------------------------------------
cat("\nU02 Salmon\n")
sh("salmon", c("index", "-t", file.path(out, "transcripts.fa"), "-i", file.path(out, "salmon_index"), "-k", "25"), file.path(out, "salmon_index.log"))
code <- run_r("later/scripts/U02_salmon_quant.R", c("--sample_sheet", file.path(out, "samples.tsv"), "--fastq_dir", file.path(out, "fastq"),
                                                    "--index", file.path(out, "salmon_index"), "--out_dir", file.path(out, "salmon"),
                                                    "--threads", "2", "--bootstraps", "0", "--jobs", "2"))
st <- if (file.exists(file.path(out, "salmon", "strandedness.txt"))) readLines(file.path(out, "salmon", "strandedness.txt")) else NA
record("Salmon quantified every sample", code == 0 && all(file.exists(file.path(out, "salmon", samples$sample, "quant.sf"))))
record("reverse-stranded library detected", identical(st, "reverse"), st)
code <- run_r("later/scripts/U02_salmon_quant.R", c("--sample_sheet", file.path(out, "samples.tsv"), "--fastq_dir", file.path(out, "nope"),
                                                    "--index", file.path(out, "salmon_index"), "--out_dir", file.path(out, "salmon_bad")))
record("missing FASTQ files stop U02", code == 1)

# ---- U03 tximport ------------------------------------------------------------------------------------------
cat("\nU03 gene counts\n")
code <- run_r("later/scripts/U03_tximport_counts.R", c("--salmon_dir", file.path(out, "salmon"), "--gtf", file.path(out, "genes.gtf"),
                                                       "--sample_sheet", file.path(out, "samples.tsv"), "--out_dir", file.path(out, "counts")))
cts <- if (file.exists(file.path(out, "counts", "counts.csv"))) read.csv(file.path(out, "counts", "counts.csv"), check.names = FALSE) else NULL
smp <- if (file.exists(file.path(out, "counts", "samples.csv"))) read.csv(file.path(out, "counts", "samples.csv")) else NULL
record("gene counts table for the pipeline", code == 0 && !is.null(cts) && setequal(cts$gene_id, c("ENSG00000000001", "ENSG00000000002")) &&
         identical(names(cts)[-1], samples$sample) && all(colSums(cts[, -1]) > 3000))
record("sample table for the pipeline", !is.null(smp) && identical(smp$sample, samples$sample) && "genotype" %in% names(smp))

# ---- U06 STAR alignment ---------------------------------------------------------------------------
cat("\nU06 STAR alignment\n")
code <- run_r("later/scripts/U06_star_align.R", c("--sample_sheet", file.path(out, "samples.tsv"), "--fastq_dir", file.path(out, "fastq"),
                                                  "--star_index", file.path(out, "star_index"), "--genome_fasta", file.path(out, "genome.fa"),
                                                  "--gtf", file.path(out, "genes.gtf"), "--read_length", "75", "--threads", "2",
                                                  "--out_dir", file.path(out, "star")))
dir.create(file.path(out, "bams"), showWarnings = FALSE)
for (s in samples$sample) for (ext in c(".bam", ".bam.bai")) file.copy(file.path(out, "star", "bams", paste0(s, ext)), file.path(out, "bams", paste0(s, ext)))
al <- if (file.exists(file.path(out, "star", "alignment_summary.csv"))) read.csv(file.path(out, "star", "alignment_summary.csv")) else NULL
record("STAR index built and every sample aligned", code == 0 && !is.null(al) && all(al$status == "success") &&
         all(file.exists(file.path(out, "bams", paste0(samples$sample, ".bam.bai")))))
record("high unique mapping rate", !is.null(al) && all(al$unique_percent > 90), if (!is.null(al)) sprintf("min %.1f%%", min(al$unique_percent)) else "")
st2 <- if (file.exists(file.path(out, "star", "strandedness.txt"))) readLines(file.path(out, "star", "strandedness.txt")) else NA
record("STAR counts also say reverse-stranded", identical(st2, "reverse"), st2)
sc <- if (file.exists(file.path(out, "star", "counts_star.csv"))) read.csv(file.path(out, "star", "counts_star.csv"), check.names = FALSE) else NULL
record("STAR gene counts for the pipeline", !is.null(sc) && setequal(sc$gene_id, c("ENSG00000000001", "ENSG00000000002")) && all(colSums(sc[, -1]) > 3000))

# ---- U04 rMATS ----------------------------------------------------------------------------------------------------
cat("\nU04 rMATS\n")
code <- run_r("later/scripts/U04_rmats.R", c("--sample_sheet", file.path(out, "samples.tsv"), "--condition_column", "genotype", "--control", "WT",
                                             "--bam_dir", file.path(out, "bams"), "--gtf", file.path(out, "genes.gtf"),
                                             "--read_length", file.path(out, "read_length", "read_length.txt"),
                                             "--strandedness", file.path(out, "salmon", "strandedness.txt"), "--read_type", "paired",
                                             "--threads", "2", "--out_dir", file.path(out, "rmats")))
cmp <- file.path(out, "rmats", "WT_vs_KO")
se <- if (file.exists(file.path(cmp, "tables", "SE.MATS.JCEC.txt"))) read.delim(file.path(cmp, "tables", "SE.MATS.JCEC.txt")) else NULL
record("rMATS ran for KO vs WT", code == 0 && !is.null(se) && nrow(se) >= 1, if (!is.null(se)) sprintf("%d SE events", nrow(se)) else "")
e2 <- if (!is.null(se)) se[se$exonStart_0base == 3000 & se$exonEnd == 3150, ] else NULL
record("E2 skipping found: less inclusion in KO", !is.null(e2) && nrow(e2) == 1 && e2$IncLevelDifference < -0.4 && e2$FDR < 0.05,
       if (!is.null(e2) && nrow(e2)) sprintf("dPSI %.2f, FDR %.2g", e2$IncLevelDifference, e2$FDR) else "event not found")
cfg <- yaml::read_yaml(file.path(cmp, "splicing_config.yaml"))
record("splicing config written", identical(cfg$group_1_name, "KO") && identical(unlist(cfg$group_2_samples), paste0("WT", 1:3)))
code <- run_r("later/scripts/U04_rmats.R", c("--sample_sheet", file.path(out, "samples.tsv"), "--condition_column", "genotype", "--control", "Control",
                                             "--bam_dir", file.path(out, "bams"), "--gtf", file.path(out, "genes.gtf"), "--read_length", "75",
                                             "--read_type", "paired", "--out_dir", file.path(out, "rmats_bad")))
record("unknown control group stops U04", code == 1)

# ---- Splicing module on the rMATS output -----------------------------------------------------------------------------
cat("\nSplicing module\n")
before <- list.dirs(out, recursive = FALSE)
code <- run_r("splicing/scripts/run_splicing.R", c("--config", file.path(cmp, "splicing_config.yaml"), "--runs_dir", out,
                                                   "--set", shQuote("run_enrichment=false")))
srun <- setdiff(list.dirs(out, recursive = FALSE), before)
sig <- if (length(srun)) read.csv(file.path(srun, "S02_splicing_results", "splicing_events_significant.csv")) else NULL
record("splicing module calls the E2 event", code == 0 && !is.null(sig) && any(sig$gene_symbol == "GENEA" & sig$event_type == "SE" & sig$direction == "down"))

# ---- U05 sashimi ---------------------------------------------------------------------------------------------------------
cat("\nU05 sashimi plots\n")
code <- run_r("later/scripts/U05_sashimi_plots.R", c("--splicing_run", srun, "--comparison_dir", cmp, "--gtf", file.path(out, "genes.gtf")))
pngs <- list.files(file.path(srun, "sashimi"), pattern = "\\.png$")
ss <- read_json_file(file.path(srun, "sashimi", "sashimi_summary.json"))
record("sashimi plot drawn with junction arcs", code == 0 && length(pngs) >= 1 && isTRUE(ss$plots$n_junctions_shown[1] >= 2),
       paste(pngs, collapse = ", "))

# ---- run_upstream.R: everything from one config ----------------------------------------------------------------------
cat("\nrun_upstream.R\n")
cfg <- list(sample_sheet = file.path(out, "samples.tsv"), fastq_dir = file.path(out, "fastq"), out_dir = file.path(out, "upstream"),
            condition_column = "genotype", control = "WT", genome_fasta = file.path(out, "genome.fa"), gtf = file.path(out, "genes.gtf"),
            star_index = file.path(out, "star_index"), salmon_index = file.path(out, "upstream_salmon_index"),
            transcripts_fasta = file.path(out, "transcripts.fa"), threads = 2, bootstraps = 0,
            steps = c("qc", "read_length", "star", "salmon", "gene_counts", "rmats", "splicing", "sashimi"),
            splicing_settings = list(run_enrichment = FALSE))
yaml::write_yaml(cfg, file.path(out, "upstream.yaml"))
code <- run_r("later/scripts/run_upstream.R", c("--config", file.path(out, "upstream.yaml")))
us <- read_json_file(file.path(out, "upstream", "upstream_summary.json"))
ok_steps <- !is.null(us) && all(vapply(us$steps, function(x) x$status, "") == "success")
record("all upstream steps succeed", code == 0 && ok_steps, if (!is.null(us)) paste(names(us$steps), collapse = ", ") else "")
usr <- us$splicing_runs[[1]]
us_sig <- if (!is.null(usr)) read.csv(file.path(usr, "S02_splicing_results", "splicing_events_significant.csv")) else NULL
record("upstream run finds the E2 event too", !is.null(us_sig) && any(us_sig$gene_symbol == "GENEA" & us_sig$direction == "down") &&
         length(list.files(file.path(usr, "sashimi"), pattern = "png$")) >= 1)
record("main-pipeline inputs listed", !is.null(us$main_pipeline_inputs$counts_salmon) && file.exists(us$main_pipeline_inputs$counts_salmon))

cat("\n================================================================\n")
n_fail <- sum(results$status == "FAIL")
cat(sprintf(" %d tests: %d passed, %d failed\n", nrow(results), sum(results$status == "PASS"), n_fail))
if (n_fail) print(results[results$status == "FAIL", ], row.names = FALSE)
cat("================================================================\n")
quit(save = "no", status = if (n_fail) 1 else 0)
