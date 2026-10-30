#!/usr/bin/env Rscript
# =============================================================================
# U06_star_align.R  --  splice-aware genome alignment with STAR (FASTQ -> BAM)
# =============================================================================
# The alignment step rMATS (U04) and the sashimi plots (U05) need. Builds the
# STAR genome index when it does not exist yet (from --genome_fasta and --gtf),
# aligns every sample, sorts and indexes the BAM files, and collects:
#   * mapping statistics per sample (STAR Log.final.out), with flags
#   * gene counts from STAR --quantMode GeneCounts (same counting rules as
#     htseq-count -m union), as a counts table for the main pipeline
#   * the library strandedness, read from those counts (independent of Salmon)
#
# Usage
#   # first time: build the index too (needs ~32 GB RAM for human / mouse)
#   Rscript later/scripts/U06_star_align.R --sample_sheet samples.tsv --fastq_dir trimmed/ \
#       --genome_fasta ref/genome.fa --gtf ref/genes.gtf --star_index ref/star_index --out_dir results/star
#   # later runs: reuse the index
#   Rscript later/scripts/U06_star_align.R --sample_sheet samples.tsv --fastq_dir trimmed/ --star_index ref/star_index --out_dir results/star
#
# Outputs (in --out_dir):
#   bams/<sample>.bam (+ .bai)       sorted, indexed alignments (input for U04 / U05)
#   star/<sample>/                   STAR's logs, splice junctions, ReadsPerGene.out.tab
#   alignment_summary.csv / .json    input reads, % uniquely mapped, % multi-mapped, % unmapped, junctions, flags
#   strandedness.txt                 unstranded | forward | reverse
#   counts_star.csv, samples.csv     gene counts + sample table -> the main pipeline
# Exit code 1 if any sample failed.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "later/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "later", "R", "later_utils.R"))
suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample_sheet", type = "character", help = "Sample sheet (required)"),
  make_option("--fastq_dir", type = "character", default = ".", help = "Folder with the FASTQ files"),
  make_option("--star_index", type = "character", help = "STAR index folder; built here if missing (required)"),
  make_option("--genome_fasta", type = "character", help = "Genome FASTA (only to build the index)"),
  make_option("--gtf", type = "character", help = "Annotation GTF (to build the index)"),
  make_option("--read_length", type = "character", default = "101", help = "Read length or U01's read_length.txt; sets sjdbOverhang [default %default]"),
  make_option("--out_dir", type = "character", help = "Output folder (required)"),
  make_option("--threads", type = "integer", default = 8L),
  make_option("--two_pass", action = "store_true", default = FALSE, help = "STAR --twopassMode Basic (better for unannotated junctions)"),
  make_option("--strandedness", type = "character", help = "unstranded | forward | reverse | a file; default: detected from the gene counts"),
  make_option("--min_unique_percent", type = "double", default = 60, help = "Flag samples with fewer uniquely mapped reads [default %default]"),
  make_option("--star", type = "character", help = "Path to STAR"),
  make_option("--samtools", type = "character", help = "Path to samtools"),
  make_option("--dry_run", action = "store_true", default = FALSE)
), description = "STAR alignment, sorted BAM files, mapping statistics and gene counts."))
for (k in c("sample_sheet", "star_index", "out_dir")) if (is_blank(opt[[k]])) fail("--", k, " is required.")
star <- need_tool("STAR", opt[["star"]])
samtools <- need_tool("samtools", opt[["samtools"]])
dir.create(file.path(opt[["out_dir"]], "bams"), recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(opt[["out_dir"]])
version <- tool_version(star)
read_len <- if (grepl("^[0-9]+$", opt[["read_length"]])) as.integer(opt[["read_length"]]) else as.integer(trimws(readLines(opt[["read_length"]], n = 1)))

# ---- Genome index ------------------------------------------------------------------------
if (!file.exists(file.path(opt[["star_index"]], "SA"))) {
  if (is_blank(opt[["genome_fasta"]]) || is_blank(opt[["gtf"]])) fail("No STAR index at ", opt[["star_index"]], ". Give --genome_fasta and --gtf to build it.")
  for (f in c(opt[["genome_fasta"]], opt[["gtf"]])) if (!file.exists(f)) fail("File not found: ", f)
  if (grepl("\\.gz$", opt[["genome_fasta"]]) || grepl("\\.gz$", opt[["gtf"]])) fail("STAR needs the genome FASTA and GTF uncompressed (gunzip them first).")
  genome_len <- file.size(opt[["genome_fasta"]])  # close enough to the genome length for the index setting
  sa_bases <- min(14L, as.integer(floor(log2(genome_len) / 2 - 1)))
  dir.create(opt[["star_index"]], recursive = TRUE, showWarnings = FALSE)
  log_info(sprintf("Building the STAR index in %s (genomeSAindexNbases %d, sjdbOverhang %d) ...", opt[["star_index"]], sa_bases, read_len - 1))
  code <- run_tool(star, c("--runMode", "genomeGenerate", "--genomeDir", shQuote(normalizePath(opt[["star_index"]])),
                           "--genomeFastaFiles", shQuote(normalizePath(opt[["genome_fasta"]])), "--sjdbGTFfile", shQuote(normalizePath(opt[["gtf"]])),
                           "--sjdbOverhang", read_len - 1, "--genomeSAindexNbases", sa_bases, "--runThreadN", opt[["threads"]],
                           "--outFileNamePrefix", shQuote(paste0(normalizePath(opt[["star_index"]]), "/"))),
                   file.path(out_dir, "logs", "star_index.log"), dry_run = opt[["dry_run"]])
  if (code != 0) fail("Building the STAR index failed; see logs/star_index.log (not enough memory is the usual cause).")
}

# ---- Align each sample ----------------------------------------------------------------------
samples <- pair_fastqs(read_sample_sheet(opt[["sample_sheet"]]), opt[["fastq_dir"]])
log_info(sprintf("%s; aligning %d samples.", version, nrow(samples)))
results <- list()
for (i in seq_len(nrow(samples))) {
  s <- samples[i, ]
  sdir <- file.path(out_dir, "star", s$sample)
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  reads <- c(shQuote(s$fastq_1), if (!is.na(s$fastq_2)) shQuote(s$fastq_2))
  gz <- grepl("\\.gz$", s$fastq_1)
  args <- c("--genomeDir", shQuote(normalizePath(opt[["star_index"]])), "--readFilesIn", reads, if (gz) c("--readFilesCommand", "gunzip", "-c"),
            "--outSAMtype", "BAM", "SortedByCoordinate", "--outSAMattributes", "NH", "HI", "AS", "nM", "XS",
            "--quantMode", "GeneCounts", "--runThreadN", opt[["threads"]], if (opt[["two_pass"]]) c("--twopassMode", "Basic"),
            "--outFileNamePrefix", shQuote(paste0(sdir, "/")))
  log_info(sprintf("%s (%d of %d) ...", s$sample, i, nrow(samples)))
  code <- run_tool(star, args, file.path(out_dir, "logs", paste0(s$sample, ".log")), dry_run = opt[["dry_run"]])
  bam_src <- file.path(sdir, "Aligned.sortedByCoord.out.bam")
  bam <- file.path(out_dir, "bams", paste0(s$sample, ".bam"))
  ok <- code == 0 && (opt[["dry_run"]] || file.exists(bam_src))
  if (ok && !opt[["dry_run"]]) {
    file.rename(bam_src, bam)
    ok <- run_tool(samtools, c("index", shQuote(bam)), file.path(out_dir, "logs", paste0(s$sample, "_index.log"))) == 0
  }
  lf <- file.path(sdir, "Log.final.out")
  stat <- function(key) {
    if (!file.exists(lf)) return(NA_real_)
    l <- grep(key, readLines(lf, warn = FALSE), fixed = TRUE, value = TRUE)
    if (!length(l)) return(NA_real_)
    suppressWarnings(as.numeric(sub("%", "", trimws(sub(".*\\|", "", l[1])))))
  }
  results[[s$sample]] <- list(sample = s$sample, status = if (ok) "success" else "failed", exit_code = code,
                              input_reads = stat("Number of input reads"), unique_percent = stat("Uniquely mapped reads %"),
                              multi_percent = stat("% of reads mapped to multiple loci"),
                              unmapped_too_short_percent = stat("% of reads unmapped: too short"),
                              splice_junctions = stat("Number of splices: Total"), bam = if (ok) bam else NA)
}
res <- do.call(rbind, lapply(results, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
res$flags <- ifelse(!is.na(res$unique_percent) & res$unique_percent < opt[["min_unique_percent"]],
                    sprintf("uniquely mapped < %g%%", opt[["min_unique_percent"]]), "")
write_table(res, file.path(out_dir, "alignment_summary.csv"))

# ---- Strandedness and gene counts --------------------------------------------------------------------
ok_samples <- res$sample[res$status == "success"]
strand <- NULL; strand_note <- "not determined"
if (length(ok_samples) && !opt[["dry_run"]]) {
  rpg <- lapply(ok_samples, function(s) {
    d <- read.delim(file.path(out_dir, "star", s, "ReadsPerGene.out.tab"), header = FALSE)
    d[!startsWith(d$V1, "N_"), ]
  })
  sums <- colSums(do.call(rbind, lapply(rpg, function(d) colSums(d[, 3:4]))))
  detected <- if (sums[1] > 4 * sums[2]) "forward" else if (sums[2] > 4 * sums[1]) "reverse" else "unstranded"
  strand <- if (!is_blank(opt[["strandedness"]])) {
    if (opt[["strandedness"]] %in% names(RMATS_LIBTYPE)) opt[["strandedness"]] else trimws(readLines(opt[["strandedness"]], n = 1))
  } else detected
  strand_note <- sprintf("Reads on the gene's strand vs the opposite strand: %s vs %s -> %s%s", format(sums[1], big.mark = ","),
                         format(sums[2], big.mark = ","), detected, if (strand != detected) sprintf(" (using the given '%s')", strand) else "")
  writeLines(strand, file.path(out_dir, "strandedness.txt"))
  col <- c(unstranded = 2, forward = 3, reverse = 4)[[strand]]
  cts <- do.call(cbind, lapply(rpg, function(d) d[[col]]))
  colnames(cts) <- ok_samples
  write_table(data.frame(gene_id = rpg[[1]]$V1, cts, check.names = FALSE), file.path(out_dir, "counts_star.csv"))
  sheet <- read_sample_sheet(opt[["sample_sheet"]])
  sheet <- sheet[!duplicated(sheet$sample) & sheet$sample %in% ok_samples, setdiff(names(sheet), c("file_name", "read_orientation")), drop = FALSE]
  write_table(sheet, file.path(out_dir, "samples.csv"))
  log_info(strand_note)
}
flagged <- res$sample[nzchar(res$flags)]
if (length(flagged)) log_warn("Low unique mapping: ", paste(flagged, collapse = ", "), ". Check species, contamination or rRNA.")
write_json_file(list(star_version = version, star_index = normalizePath(opt[["star_index"]]), two_pass = opt[["two_pass"]],
                     strandedness = strand, strandedness_note = strand_note, flagged = as.list(flagged),
                     next_steps = c("U04_rmats.R --bam_dir <out_dir>/bams --strandedness <out_dir>/strandedness.txt",
                                    "Main pipeline: counts_star.csv + samples.csv"),
                     samples = unname(results)), file.path(out_dir, "alignment_summary.json"))
n_fail <- sum(res$status == "failed")
if (n_fail) fail(sprintf("Alignment failed for %d sample(s): %s", n_fail, paste(res$sample[res$status == "failed"], collapse = ", ")))
log_info(sprintf("Aligned %d samples (median %.1f%% uniquely mapped) -> %s", nrow(res), median(res$unique_percent, na.rm = TRUE), out_dir))
