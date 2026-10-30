# later/: from FASTQ files to counts, splicing and isoforms (not used by the website yet)

The website currently starts from a **counts table**, so these steps are not needed yet. They are kept here, tested, for when raw sequencing files arrive through the **Globus upload**. They include R versions of the lab's Python scripts (`read_length.py`, `salmon_alignment.py`, `rmats_analysis.py`, `rMATS_ggsashimi.py`), with the problems listed at the end fixed, plus the steps those scripts assumed but did not include: read QC and trimming, and STAR alignment.

```
FASTQ ──U00 fastp QC + trimming──> trimmed FASTQ
          │
          ├──U01 read length
          │
          ├──U06 STAR alignment ──> BAM ──U04 rMATS ──> splicing module (splicing/) ──U05 sashimi plots
          │        └──> STAR gene counts + strandedness ─────────────┐
          │                                                          ├──> main pipeline (scripts/run_pipeline.R)
          └──U02 Salmon ──U03 tximport──> counts.csv + samples.csv ──┘
                   └──> transcript quantification ──> isoform module (isoforms/): DEI + DUI
```

**One command for everything:** `Rscript later/scripts/run_upstream.R --config upstream.yaml`. It runs the steps you list in order, stops at the first failure, and writes `upstream_summary.json` with every output path. This is what the Globus integration will call.

## Scripts

| Script | From | Does | Writes |
|---|---|---|---|
| `U00_fastq_qc_trim.R` | new | fastp: adapter detection and removal, quality trimming of read ends, removal of too-short reads, a QC report per sample; flags samples that keep few reads or have few Q30 bases | trimmed FASTQ, `reports/<sample>.html`, `qc_summary.csv`, `samples_trimmed.tsv` |
| `U01_read_length.R` | `read_length.py` | Read length from 10,000 reads per file (not only the first read); flags variable lengths after trimming | `read_length.txt`, `read_length.json` |
| `U06_star_align.R` | new | STAR splice-aware alignment; builds the genome index the first time; sorted and indexed BAM files; mapping statistics with flags; gene counts (`--quantMode GeneCounts`); detects strandedness from those counts | `bams/*.bam`, `alignment_summary.csv`, `counts_star.csv`, `samples.csv`, `strandedness.txt` |
| `U02_salmon_quant.R` | `salmon_alignment.py` | Salmon quantification with bias corrections and bootstraps; paired-end and single-end; samples in parallel; checks every exit code; detects strandedness; flags low mapping rates | `<sample>/quant.sf`, `salmon_run_summary.csv/json`, `strandedness.txt` |
| `U03_tximport_counts.R` | new (bridge) | tximport gene counts (`lengthScaledTPM`) with the transcript-to-gene mapping from the GTF | `counts.csv` and `samples.csv`: the main pipeline's two input files |
| `U04_rmats.R` | `rmats_analysis.py` | rMATS for each group against the control; valid rMATS library types; variable read length from U01 | `<control>_vs_<group>/tables/*.MATS.JCEC.txt`, `splicing_config.yaml` |
| `U05_sashimi_plots.R` | `rMATS_ggsashimi.py` + plotting part of `rmats_analysis.py` | Sashimi plots (mean coverage, junction arcs with read counts, transcripts, shaded event exon) for the top significant events of a splicing run | `<splicing run>/sashimi/*.png` |
| `run_upstream.R` | new | Runs all of the above, plus the splicing and isoform modules, from one config | `upstream_summary.json` |

The lab's `rMATS_Plots.R` (summary figures) and GSE script (pathways of spliced genes) are covered by the [splicing module](../splicing/README.md). Its S02 step draws event counts, the ΔPSI volcano, the ΔPSI distribution, the PSI scatter, the per-sample PSI heatmap, and splicing vs expression. Its S03 step tests pathways for all spliced genes, for genes that only gain or only lose inclusion, and per event type. **DEI** (differential isoform expression) and **DUI** (differential isoform usage) are in the [isoform module](../isoforms/README.md).

Each script prints its usage with `--help`.

## One config for everything

```yaml
# upstream.yaml
sample_sheet: samples.tsv
fastq_dir: fastq/
out_dir: results/
condition_column: genotype
control: WT
genome_fasta: ref/GRCm39.primary_assembly.genome.fa   # uncompressed; STAR index built once
gtf: ref/gencode.vM33.annotation.gtf
star_index: ref/star_index
salmon_index: ref/salmon_index                         # built from transcripts_fasta when missing
transcripts_fasta: ref/gencode.vM33.transcripts.fa
threads: 16
steps: [qc, read_length, star, salmon, gene_counts, rmats, splicing, sashimi, isoforms]
splicing_settings: {delta_psi_cutoff: 0.1}
isoform_settings: {min_delta_usage: 0.1}
```

```bash
Rscript later/scripts/run_upstream.R --config upstream.yaml
```

The sample sheet has one row per FASTQ file, in TSV or CSV format; single-end samples have one row with `R1` (or `SE`):

```
file_name          sample  genotype  read_orientation  batch
WT1_R1.fastq.gz    WT1     WT        R1                1
WT1_R2.fastq.gz    WT1     WT        R2                1
```

## Requirements

- **Tools:** fastp, STAR, samtools, Salmon and rMATS 4.x. One conda environment covers them all: `conda create -n rnaseq-later -c conda-forge -c bioconda fastp star samtools salmon rmats`.
- **R packages:** the code bank's packages plus `tximport`, and `DEXSeq` and `stageR` for isoform usage (all on Bioconductor). U05 needs only ggplot2; it reads BAM files through samtools.
- **Reference files:** genome FASTA, GTF and transcript FASTA from the same release (e.g. GENCODE human v44 / GRCh38, mouse M33 / GRCm39).
- **Hardware:** FASTQ and BAM files are tens of GB, and building or loading a human or mouse STAR index needs about 32 GB of memory. These steps belong on the server (Ada), not a laptop.

## What changed from the Python versions

- `rmats_gse_preprocessor.py` filled in random fold changes and p-values (`np.random`) for spliced genes before the enrichment step. Any result from those numbers is meaningless and changes every run. Replaced by over-representation of the gene lists (splicing module, S03).
- The sashimi window was taken from the minimum and maximum of every numeric column of the rMATS table, including ID and FDR, so most plots showed a window starting near 0. U05 uses the event's own exon coordinates.
- Almost every step was wrapped in `try/except: pass`, rMATS's exit code was ignored, and the script exited 0 on errors. Now every tool's exit code is checked, and a failure stops the script with exit code 1 and a message.
- The unstranded fallback was `unstranded`, which rMATS does not accept. It is now `fr-unstranded`. Strandedness is now detected twice, from Salmon and from STAR's counts, instead of coming from a file made elsewhere.
- When nothing was significant, the first annotated events were plotted as if they were results. Now nothing is plotted, and the summary says why.
- Salmon accepted only paired-end samples; single-end samples now work too.
- The read length came from the first read only; U01 checks 10,000 reads per file.

## Test

```bash
Rscript later/tests/run_later_tests.R     # needs fastp, salmon, STAR, samtools, rmats.py on the PATH and tximport
```

It builds a small synthetic genome with a cassette exon that WT includes in 80% of transcripts and KO in 20%, simulates reverse-stranded paired-end reads, and runs every script, then everything again through `run_upstream.R`. It checks that:

- fastp keeps all reads;
- the read length is 75;
- Salmon and STAR both detect the reverse-stranded library;
- gene counts come out from both Salmon and STAR;
- rMATS finds the exon skipping (ΔPSI ≈ −0.58);
- the splicing module calls the event;
- the sashimi plot is drawn.

There are 22 tests and they take about a minute. Without the tools installed, it says which are missing and exits.
