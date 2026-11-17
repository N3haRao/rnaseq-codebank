#!/usr/bin/env Rscript
# =============================================================================
# majiq/tests/run_majiq_tests.R  --  tests of the MAJIQ module
# =============================================================================
# Run from the repository root:  Rscript majiq/tests/run_majiq_tests.R
# Output: majiq/tests/output/.
#
# MAJIQ needs an academic licence and real BAMs, so these tests use stand-in
# majiq / majiq-v3 / voila programs (shell scripts written below) that record their
# command lines and write voila-style tables. They check that:
#   - the commands match the lab's MAJIQ v3 workflow (build, psi-coverage, sg-coverage,
#     deltapsi with the reference as -psi1, voila tsv / modulize with the thresholds);
#   - the results step reads voila's tables (JSON comment block, ';'-separated junction
#     values), finds the changing LSVs with the right direction, and counts modulize events;
#   - missing MAJIQ, licence, BAMs, BAM index or a GTF instead of GFF3 fail with clear checks.
# With MAJIQ installed (MAJIQ_TEST_REAL=1), the real programs are only version-checked.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "majiq/tests/x"), "..", ".."))
})
setwd(REPO_ROOT)
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "majiq", "R", "majiq_utils.R"))

out <- file.path(REPO_ROOT, "majiq", "tests", "output")
unlink(out, recursive = TRUE)
dir.create(out, recursive = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")
results <- data.frame(test = character(0), status = character(0), detail = character(0))
record <- function(test, ok, detail = "") {
  if (!length(detail) || is.na(detail[1])) detail <- "" else detail <- as.character(detail[1])
  results[nrow(results) + 1, ] <<- list(test, if (isTRUE(ok)) "PASS" else "FAIL", detail)
  cat(sprintf("  [%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", test, if (nzchar(detail)) paste0("- ", detail) else ""))
}

# ---- Stand-in MAJIQ ----------------------------------------------------------------------------
stub <- file.path(out, "fake_majiq_bin"); dir.create(stub)
calls <- file.path(out, "calls.log")
writeLines(c("#!/usr/bin/env bash", sprintf('echo "majiq $*" >> "%s"', calls),
  'arg() { local k=$1; shift; while [ $# -gt 0 ]; do [ "$1" = "$k" ] && { echo "$2"; return; }; shift; done; }',
  'case "$1" in',
  '  --version) echo "majiq-v3 v3.0.23";;',
  '  build) shift; while [ "$1" = "--license" ]; do shift 2; done; gff=$1; exp=$2; outd=$3; mkdir -p "$outd/splicegraph.zarr";',
  '         tail -n +2 "$exp" | while IFS=$(printf "\\t") read -r p g; do b=$(basename "$p" .bam); echo sj > "$outd/$b.sj"; done;;',
  '  psi-coverage) shift; while [ "$1" = "--license" ]; do shift 2; done; echo psicov > "$2";;',
  '  deltapsi) echo dpsi > "$(arg --output-voila "$@")"; echo "# {}" > "$(arg --output-tsv "$@")";;',
  '  *) echo "unknown $1" >&2; exit 2;;',
  'esac'), file.path(stub, "majiq"))
writeLines(c("#!/usr/bin/env bash", sprintf('echo "majiq-v3 $*" >> "%s"', calls),
  '[ "$1" = sg-coverage ] && { shift; while [ "$1" = "--license" ]; do shift 2; done; echo sgc > "$2"; exit 0; }; exit 2'),
  file.path(stub, "majiq-v3"))
# voila: 100 LSVs; LSV 1-6 more inclusion in KO (dPSI +0.4), 7-10 less (-0.35), all with P 0.99; the rest P 0.05
writeLines(c("#!/usr/bin/env bash", sprintf('echo "voila $*" >> "%s"', calls),
  'arg() { local k=$1; shift; while [ $# -gt 0 ]; do [ "$1" = "$k" ] && { echo "$2"; return; }; shift; done; }',
  'mode=$3',
  'if [ "$mode" = tsv ]; then f=$(arg -f "$@"); changing_only=0; for a in "$@"; do [ "$a" = --probability-threshold ] && changing_only=1; done',
  '  { echo "# {"; echo "#   \\"voila_version\\": \\"3.0.23\\""; echo "# }"',
  '    if [ "${STUB_COLUMNS:-}" = odd ]; then printf "gene\\tthing\\tvalue\\n"; printf "G\\tx\\t1\\n"; else',
  '    printf "gene_name\\tgene_id\\tlsv_id\\tmean_dpsi_per_lsv_junction\\tprobability_changing\\tprobability_non_changing\\tWT_mean_psi\\tKO_mean_psi\\tlsv_type\\tnum_junctions\\tnum_exons\\tde_novo_junctions\\tseqid\\tstrand\\tjunctions_coords\\n"',
  '    for i in $(seq 1 100); do',
  '      if [ $i -le 6 ]; then d="0.400;-0.380"; p="0.990;0.985"; elif [ $i -le 10 ]; then d="-0.350;0.340"; p="0.992;0.980"; else d="0.020;-0.010"; p="0.050;0.040"; fi',
  '      [ $changing_only = 1 ] && [ $i -gt 10 ] && continue',
  '      printf "GENE%03d\\tENSG%011d.1\\tENSG%011d.1:s:%d-%d\\t%s\\t%s\\t0.95;0.96\\t0.40;0.60\\t0.80;0.20\\ts|1e1.1o1|1e2.1o1\\t2\\t3\\t0\\tchr1\\t+\\t%d-%d;%d-%d\\n" $i $i $i $((i*1000)) $((i*1000+100)) "$d" "$p" $((i*1000+100)) $((i*1000+400)) $((i*1000+100)) $((i*1000+800))',
  '    done; fi; } > "$f"',
  'elif [ "$mode" = modulize ]; then d=$(arg -d "$@"); mkdir -p "$d"',
  '  { echo "# {}"; printf "module_id\\tgene_id\\tevent_changing\\n"; printf "m1\\tg1\\tTrue\\nm2\\tg2\\tFalse\\nm3\\tg3\\tTrue\\n"; } > "$d/cassette.tsv"',
  '  { echo "# {}"; printf "module_id\\tgene_id\\tevent_changing\\n"; printf "m4\\tg4\\tTrue\\nm5\\tg5\\tFalse\\n"; } > "$d/alt3prime.tsv"',
  '  { echo "# {}"; printf "module_id\\tgene_id\\n"; printf "m1\\tg1\\n"; } > "$d/summary.tsv"',
  'else exit 2; fi'), file.path(stub, "voila"))
Sys.chmod(file.path(stub, c("majiq", "majiq-v3", "voila")), "0755")
writeLines("fake academic licence for tests", file.path(out, "majiq_license_test.lic"))

# ---- Inputs: three BAM layouts, a GFF3 ------------------------------------------------------------
samples <- c(paste0("WT", 1:3), paste0("KO", 1:3))
bd <- file.path(out, "bams"); dir.create(bd)
for (s in samples[1:2]) { writeLines("bam", file.path(bd, paste0(s, ".bam"))); writeLines("bai", file.path(bd, paste0(s, ".bam.bai"))) }
for (s in samples[3:4]) { writeLines("bam", file.path(bd, paste0(s, ".sorted.bam"))); writeLines("bai", file.path(bd, paste0(s, ".sorted.bam.bai"))) }
for (s in samples[5:6]) { dir.create(file.path(bd, s)); writeLines("bam", file.path(bd, s, "Aligned.sortedByCoord.out.bam"))
                          writeLines("bai", file.path(bd, s, "Aligned.sortedByCoord.out.bam.bai")) }
writeLines(c("##gff-version 3", "chr1\tHAVANA\tgene\t1000\t9000\t.\t+\t.\tID=ENSG00000000001.1;gene_id=ENSG00000000001.1;gene_name=GENE001",
             "chr1\tHAVANA\ttranscript\t1000\t9000\t.\t+\t.\tID=ENST00000000001.1;Parent=ENSG00000000001.1"), file.path(out, "genes.gff3"))
write.csv(data.frame(sample = samples, genotype = rep(c("WT", "KO"), each = 3)), file.path(out, "samples.csv"), row.names = FALSE)

run_mj <- function(args, env = character(0)) {
  before <- list.dirs(out, recursive = FALSE)
  log <- tempfile(fileext = ".log")
  code <- system2(rscript, c("majiq/scripts/run_majiq.R", "--runs_dir", shQuote(out), args), stdout = log, stderr = log,
                  env = c("MAJIQ_BIN_DIR=", "MAJIQ_LICENSE_FILE=", "COPILOT_HOME=", env))
  if (code != 0) cat(paste0("    | ", tail(readLines(log), 8), "\n"), sep = "")
  list(code = code, run_dir = setdiff(list.dirs(out, recursive = FALSE), before)[1])
}
common <- c("--bam_dir", shQuote(bd), "--metadata_file", shQuote(file.path(out, "samples.csv")), "--gff3_file", shQuote(file.path(out, "genes.gff3")),
            "--condition", "genotype", "--reference_level", "WT", "--test_level", "KO",
            "--set", shQuote(sprintf("majiq_bin_dir=%s;majiq_license_file=%s;threads=4", stub, file.path(out, "majiq_license_test.lic"))))

# ---- 1. Full run with the stand-in MAJIQ ------------------------------------------------------------
cat("Full run (stand-in MAJIQ):\n")
r <- run_mj(common)
summ <- read_json_file(file.path(r$run_dir, "run_summary.json"))
record("run finishes with status success", r$code == 0 && identical(summ$status, "success"), summ$status %||% "")
cl <- readLines(calls)
has <- function(pattern) any(grepl(pattern, cl, fixed = TRUE))
record("BAMs of all three layouts linked as <sample>.bam",
       all(file.exists(file.path(r$run_dir, "MJ01_validate_majiq", "bams", paste0(samples, c(".bam")))) ) &&
         all(file.exists(file.path(r$run_dir, "MJ01_validate_majiq", "bams", paste0(samples, ".bam.bai")))))
exp <- read.delim(file.path(r$run_dir, "MJ01_validate_majiq", "experiments.tsv"))
record("experiments.tsv: path + group", identical(names(exp), c("path", "group")) && identical(exp$group, rep(c("WT", "KO"), each = 3)))
record("build: licence, GFF3, strandness AUTO, threads", has("majiq build --license") && has("--strandness AUTO --nthreads 4"))
record("psi-coverage per group from its 3 junction files",
       sum(grepl("majiq psi-coverage", cl)) == 2 && has("psicov/WT.psicov") && has("WT1.sj") && has("KO3.sj"))
record("sg-coverage with majiq-v3 (for modulize)", sum(grepl("majiq-v3 sg-coverage", cl)) == 2)
record("deltapsi: reference first (-psi1 WT), -n WT KO, changing threshold 0.2",
       any(grepl("majiq deltapsi .*-psi1 \\S*WT\\.psicov -psi2 \\S*KO\\.psicov -n WT KO .*--changing-threshold 0\\.2", cl)))
record("voila tsv: all LSVs with read counts", has("--show-all --show-read-counts"))
record("voila tsv: changing LSVs (0.2, P 0.95)", has("--changing-between-group-dpsi 0.2 --probability-threshold 0.95"))
record("voila modulize: sgc folder, thresholds, decomplexify 0.05",
       any(grepl("voila --license \\S+ modulize .* \\S*/sgc -d .*--probability-changing-threshold 0\\.95 --decomplexify-psi-threshold 0\\.05", cl)))
ms <- read_json_file(file.path(r$run_dir, "MJ05_majiq_results", "majiq_summary.json"))
record("10 changing LSVs found (6 more, 4 less inclusion)", identical(as.integer(ms$n_lsvs_changing), 10L) && ms$n_up == 6 && ms$n_down == 4,
       sprintf("%s (%s up, %s down)", ms$n_lsvs_changing, ms$n_up, ms$n_down))
record("agrees with voila's own changing table", identical(as.integer(ms$n_lsvs_changing_voila), 10L))
lsv <- data.table::fread(file.path(r$run_dir, "MJ05_majiq_results", "majiq_lsvs_changing.csv"), data.table = FALSE)
record("junction with the highest probability is reported", all(abs(lsv$dpsi) >= 0.34) && all(lsv$probability_changing >= 0.98))
record("modulize events counted per type", identical(ms$modulize, "done") && any(ms$modulize_by_type$event_type == "cassette" & ms$modulize_by_type$n_events == 3 &
                                                                                ms$modulize_by_type$n_changing == 2))
for (f in c("lsv_volcano.png", "event_types.png")) record(paste("figure", f), file.exists(file.path(r$run_dir, "MJ05_majiq_results", f)))
record("methods text written", any(grepl("MAJIQ 3.0.23", readLines(file.path(r$run_dir, "MJ06_majiq_report", "methods.md")))))

# ---- 2. Without modulize ----------------------------------------------------------------------------
writeLines(character(0), calls)
r2 <- run_mj(c(common[-length(common)], shQuote(sprintf("majiq_bin_dir=%s;majiq_license_file=%s;run_modulize=false", stub, file.path(out, "majiq_license_test.lic")))))
record("run_modulize false: no sg-coverage, no modulize", r2$code == 0 && !any(grepl("sg-coverage|modulize", readLines(calls))))

# ---- 3. voila columns that are not recognised -----------------------------------------------------------
r3 <- run_mj(common, env = "STUB_COLUMNS=odd")
s3 <- read_json_file(file.path(r3$run_dir, "run_summary.json"))
log5 <- readLines(file.path(r3$run_dir, "logs", "MJ05_majiq_results.log"))
record("unknown voila columns: results step stops and lists the columns", identical(s3$status, "success_with_errors") &&
         any(grepl("Columns in the file: gene, thing, value", log5)))

# ---- 4. Bad inputs -----------------------------------------------------------------------------------
cat("\nBad inputs:\n")
check_fail <- function(run, id) {
  rep <- read_json_file(file.path(run$run_dir, "MJ01_validate_majiq", "validation_report.json"))
  run$code == 1 && id %in% rep$checks$id[rep$checks$status == "fail"]
}
nob <- c(common[1:(length(common) - 1)], shQuote(sprintf("majiq_bin_dir=%s;majiq_license_file=%s", file.path(out, "nowhere"), file.path(out, "majiq_license_test.lic"))))
record("missing MAJIQ fails with a clear message", check_fail(run_mj(nob, env = "PATH=/usr/bin:/bin"), "majiq"))
nol <- c(common[1:(length(common) - 1)], shQuote(sprintf("majiq_bin_dir=%s;majiq_license_file=%s", stub, file.path(out, "none.lic"))))
lic_dir <- file.path(out, "nolic_home"); dir.create(lic_dir)
record("missing licence fails", check_fail(run_mj(nol, env = paste0("HOME=", lic_dir)), "majiq_license"))
invisible(file.remove(file.path(bd, "KO3", "Aligned.sortedByCoord.out.bam.bai")))
record("BAM without index fails (no samtools here)", if (nzchar(Sys.which("samtools"))) TRUE else check_fail(run_mj(common, env = "PATH=/usr/bin:/bin"), "bam_index"))
writeLines("bai", file.path(bd, "KO3", "Aligned.sortedByCoord.out.bam.bai"))
writeLines('chr1\tT\texon\t1\t100\t.\t+\t.\tgene_id "G"; transcript_id "T";', file.path(out, "genes.gtf"))
record("a GTF instead of GFF3 fails", check_fail(run_mj(sub("genes.gff3", "genes.gtf", common, fixed = TRUE)), "gff3_file"))
unlink(file.path(bd, "KO2"), recursive = TRUE)
record("a sample without BAM fails", check_fail(run_mj(common), "bam_files"))

if (identical(Sys.getenv("MAJIQ_TEST_REAL"), "1")) {
  bin <- find_majiq_bin(list())
  record("real MAJIQ found and starts", !is.null(bin) && system2(file.path(bin, "majiq"), "--version", stdout = FALSE, stderr = FALSE) == 0, bin %||% "not found")
}

cat("\n================================================================\n")
n_fail <- sum(results$status == "FAIL")
cat(sprintf(" %d tests: %d passed, %d failed\n", nrow(results), sum(results$status == "PASS"), n_fail))
if (n_fail) print(results[results$status == "FAIL", ], row.names = FALSE)
cat("================================================================\n")
quit(save = "no", status = if (n_fail) 1 else 0)
