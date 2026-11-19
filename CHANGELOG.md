# Changelog

Each release is tagged in git (`git tag v0.1.0`). Every run records the commit it used in `provenance.json`.

## v0.7.2

- **SUPPA2 (SU03): the reference group's PSI was missing.** SUPPA2 does not name the columns of its `.psivec` file after the samples but `<psi file name up to the first dot>_1, _2, ...` (e.g. `Non.targeting.psi` gives `Non_1`, `Non_2`, ...). SU03 matched those columns to sample IDs, so a group's PSI was found only when its samples happened to be named `<group>_1`, `<group>_2`, ...; otherwise `psi_reference` / `psi_test` were empty and `psi_scatter.png` was not made. dPSI and p-values (SUPPA's own) were not affected. SU03 now names the columns by position from the `.psi` headers (reference samples, then test samples) and stops if the counts differ.
- `suppa/tests/run_suppa_tests.R`: the synthetic samples are now named unlike their groups (`mouse1-wt`, ...), and a new test checks both groups' PSI are reported and dPSI = PSI(test) − PSI(reference) (20 tests).

## v0.7.1

- `scripts/release.sh` picks an R that can run the tests: the server's `Rscript` when it has the code bank's packages, otherwise the site image's (new `scripts/rscript_image.sh`, which runs the image's Rscript with `$COPILOT_HOME` bound), otherwise it stops with a clear message. On fuji the system R lacks optparse, fgsea and msigdbr, so every module test failed at start-up. `RSCRIPT=...` still overrides.

## v0.7.0

- New **SUPPA2 module** (`suppa/`): alternative splicing from Salmon transcript abundance, adapted from the lab's SUPPA2 scripts. Built like the other modules (run folders, logs, provenance, manifest entry with a `site` block, tests).
  - SU01 checks SUPPA2, the groups, quant.sf per sample and the GTF match, and writes SUPPA's TPM table per group. SU02 runs `generateEvents` (SE, SS, MX, RI, FL; strict) and `psiPerEvent`. SU03 runs `diffSplice -m empirical` with the reference first (dPSI = test − reference), gene-level correction by default and `-pa` only for paired samples, and writes event and gene tables, a summary and figures. SU04 (optional, `rmats_dir`) matches events with rMATS by their splice junctions, ported from the lab's `compare_SUPPA_rMATS.py` (MX sign on the − strand handled). SU05 writes the methods.
  - `suppa/install/install_suppa.sh`: SUPPA2 in a conda environment at `$COPILOT_HOME/envs/suppa`, with pandas < 2 and statsmodels 0.13 pinned (newer versions break diffSplice).
  - `suppa/tests/run_suppa_tests.R`: 19 tests with real SUPPA2 on a synthetic annotation and Salmon output with planted skipped exons and retained introns.
- New **MAJIQ module** (`majiq/`): local splicing variations from BAMs with MAJIQ v3, adapted from the lab's `run_MAJIQ.sh`: build, psi-coverage per group, sg-coverage, deltapsi (reference as `-psi1`), voila tsv (all and changing) and voila modulize; MJ05 turns voila's tables into changing LSVs, genes, event types, a summary and figures.
  - MAJIQ has an academic licence and is not in the site image: `majiq/install/install_majiq.sh` installs it per server at `$COPILOT_HOME/envs/majiq` with the licence file in `$COPILOT_HOME` (no personal details in the script).
  - `majiq/tests/run_majiq_tests.R`: 24 tests with stand-in majiq / voila programs (commands, table parsing, failures); `MAJIQ_TEST_REAL=1` also checks a real installation.
  - voila's column names vary between versions: MJ05 finds them by name and, if one is missing, stops with the columns it saw (please send that list).
- Methods policy: SUPPA2 and MAJIQ moved from `not_allowed` to `allowed` (`alternative_splicing_suppa`, `alternative_splicing_majiq`); leafcutter, JunctionSeq, VAST-tools, SplAdder, Whippet and SGSeq stay not allowed.
- Manifest: new optional keys `install` (install script) and `licence` per module.
- `scripts/release.sh` tests every module listed in the manifest that changed (it had a fixed list without the new modules).

### For the copilot (site changes this release needs; the release agent's job)

The `site` blocks use only existing handlers (`sample_folders`, `file`, `files`), so both modules appear on the "More analyses" tab with upload forms as they are. On a server with "Start from FASTQ", runs made from FASTQ results should also be able to use those results directly, as `splicing` and `isoforms` already do (`backend/app/fastq.py`, `extras_inputs` / `extras_available` / `_auto_fastq_extras`):
- **suppa:** `quant_zip` = the job's `salmon/<sample>/quant.sf` (as for isoforms); `gtf` = the reference GTF the job used (`refs/<species>/...annotation.gtf`); `rmats_tables` = the job's rMATS tables for this comparison when present (as for splicing). Cheap enough to start automatically after the analysis, like isoforms.
- **majiq:** `bams` = the job's STAR BAMs (`star/star/<sample>/Aligned.sortedByCoord.out.bam` + `.bai`; link them rather than copy, they are tens of GB) and `gff3` = a GFF3 of the job's reference release. The references (`REFERENCES` in fastq.py) need a GFF3 entry per species (GENCODE: `gencode.v44.primary_assembly.annotation.gff3.gz`, `gencode.vM33.primary_assembly.annotation.gff3.gz`; MAJIQ reads .gz). Only offer MAJIQ when it is installed (`$COPILOT_HOME/envs/majiq/bin/majiq` and a licence file); do not start it automatically (hours, licence).
- No change to the site image: SUPPA2 and MAJIQ live in `$COPILOT_HOME/envs/` (bound into the container); `COPILOT_HOME` is already set inside it.

## v0.6.0

- Every module's manifest entry now has a **`site` block** that describes how the website offers it: uploads (with a named handler), settings the user can change, values taken from the run, result tiles, headline sentences, figures, top tables, downloads, methods and the files the chat assistant reads. The copilot draws its "More analyses" tab, report section, PowerPoint slide, chat context and code .zip from it, so a new module, setting, figure or table needs no website code. `compare` is marked as its own page, `later` as not on the site yet.
- New `docs/ADDING_A_MODULE.md`: how to add a module and every `site` field.
- New `tests/check_manifest.R`: checks that every module's parameters match its defaults file, that its scripts exist, and that its `site` block only uses known building blocks. Run it before every release.
- Manifest: added the parameters that were in the defaults but not listed (`run_name` and `seed` for splicing, isoforms and compare; `organism`, `n_label` and `group_colours` for isoforms).
- Splicing S02: `splicing_vs_expression.png` is no longer drawn (empty) when no gene matches the DE results; a warning says why.

## v0.5.0

- New **compare module** (`compare/`): compares the results of 2 to 8 RNA-seq studies, never their raw counts. Each study keeps its own normalisation and design. Built like the other modules (run folders with copies of the inputs, parameters, logs, provenance, manifest entry, tests).
  - Studies are main-pipeline runs or published DE tables (DESeq2, edgeR or limma column names are recognised); `flip: true` reverses a study that compared the groups the other way round.
  - C01 matches genes by gene ID, Entrez ID or symbol, or human-mouse one-to-one orthologs (babelgene), and calls significant genes with one rule for every study.
  - C02: fold-change correlation (Pearson, Spearman), overlap of significant genes by direction with one-sided Fisher's exact tests, Jaccard index, direction agreement, scatter and overlap figures, correlation heatmap, genes significant in several studies, intersections (3+ studies). Warns when fold changes are negatively correlated (a study probably needs `flip`).
  - C03: stratified rank-rank hypergeometric overlap (RRHO2) with Benjamini-Yekutieli correction, maps and the genes at the strongest cells.
  - C04: pathway comparison from each run's own GSEA (GSEA run here for tables): NES heatmap across studies, NES correlation and shared pathways per pair.
  - C05: methods paragraph.
  - `compare/tests/run_compare_tests.R`: 32 tests with planted concordant, study-specific and discordant genes, a flipped edgeR table, and a human vs mouse comparison.
- Methods policy: cross-study comparison added as allowed; merging counts across studies and other RRHO implementations listed under `not_allowed`.
- `babelgene` added to `scripts/00_setup_environment.R`.

## v0.4.0 

- New **isoform module** (`isoforms/`), built like the other modules (run folders, logs, provenance, manifest entry, tests):
  - DEI, differential isoform expression: DESeq2 per transcript, with the main pipeline's filter, design, test and shrinkage, on tximport counts with length offsets.
  - DUI, differential isoform usage: DEXSeq on `dtuScaledTPM` transcript counts within genes, with a DRIMSeq-style filter and stageR two-stage FDR control (Love, Soneson & Patro 2018). It reports isoform switches and draws transcript-share figures.
  - `isoforms/tests/run_isoform_tests.R`: 15 tests with planted expression-only and switch genes.
- `later/` additions:
  - `U00_fastq_qc_trim.R`: fastp trimming and QC, plus a sample sheet for the trimmed files.
  - `U06_star_align.R`: STAR index build and alignment, sorted and indexed BAM files, mapping statistics, gene counts, and strandedness detected from them.
  - `run_upstream.R`: every step, including the splicing and isoform modules, from one config, with a summary of outputs.
  - The test now covers all of them (22 tests).
- Splicing module, covering the lab's `rMATS_Plots.R` and GSE script:
  - a ΔPSI distribution figure;
  - optional `de_run`, which sorts genes into splicing only / expression only / both, with a figure;
  - S03 also tests genes that only gain or only lose inclusion.
  - 20 tests.
- `R/gtf.R`: GTF reading and transcript-to-gene mapping shared by `later/` and `isoforms/`.
- Methods policy: added isoform expression (DESeq2), isoform usage (DEXSeq + stageR) and read processing (fastp, STAR, Salmon) as allowed. DEXSeq moved from not allowed to allowed (for DUI). Other isoform-usage tools (DRIMSeq, IsoformSwitchAnalyzeR, satuRn, swish) listed under `not_allowed`.
- All `later/` scripts read their options exactly (`opt[["x"]]`). R's `$` partial matching had made `opt$star` return `--star_index`.

## v0.3.0 

- New **splicing module** (`splicing/`): differential alternative splicing from rMATS output tables, built like the main pipeline (run folders, parameters, logs, provenance, run summary, manifest entry under `modules`, tests).
  - S01 checks the tables (columns, statistics present, replicates, group names) and standardises all five event types (1-based coordinates).
  - S02 calls significant events (rMATS FDR, |dPSI| and a read-coverage filter in both groups), writes tables and figures (event counts, dPSI volcano, PSI scatter, per-sample PSI heatmap).
  - S03 tests pathway over-representation among spliced genes with `fgsea::fora` against all genes with tested events (no GSEA: splicing has no gene-level ranking).
  - S04 writes a methods paragraph.
  - `splicing/tests/run_splicing_tests.R`: 16 tests on synthetic rMATS tables with planted events.
- New **`later/` folder**: R versions of the lab's FASTQ / BAM scripts, kept for the Globus upload integration and not used by the website yet.
  - U01 read length (checks many reads, flags variable lengths).
  - U02 Salmon (paired and single-end, checks every exit code, detects strandedness).
  - U03 tximport gene counts plus a sample table for the main pipeline.
  - U04 rMATS for each group against the control, with valid library types; it writes a config for the splicing module.
  - U05 sashimi plots (samtools + ggplot2, windows from the event's own coordinates).
  - `later/tests/run_later_tests.R` runs all of them end to end on a synthetic genome (needs salmon, STAR, samtools, rMATS).
- Fixed from the Python versions: enrichment input no longer uses random fold changes and p-values; sashimi windows come from the event's own coordinates; failures stop with an error instead of exiting 0; the unstranded library type is `fr-unstranded`.
- `R/genesets.R`: gene-set loading moved out of step 06 so both enrichment steps share it (step 06 output unchanged).
- Methods policy: rMATS added as the splicing method; other splicing tools (leafcutter, MAJIQ, SUPPA2, DEXSeq, ...) listed under `not_allowed`.

## v0.2.0 

- Plot options (they change how figures look, never what is computed): `plot_palette` (okabe_ito, set2, viridis, greys), `group_colours` (e.g. `["WT=#000000"]`), `up_colour` / `down_colour` for volcano and MA plots, `volcano_label_genes` (label exactly these genes), `pca_colour_by`, `pca_shape_by` (`none` for no shapes) and `pca_label_samples`. The colours used are recorded in `plots_summary.json` and `exploratory_summary.json`.
- New `methods` section in `manifest/tools.yaml`: the statistical methods of the code bank (`allowed`) and common alternatives agents must not substitute (`not_allowed`, with what to use instead). The website's guardrails read it.
- New parameter types in the manifest: `column` (a sample-table column) and `colour`.
- Test H: plot options change the figures but leave the DE table identical.

## Earlier changes since v0.1.0

- Added `scripts/07_report.R`, listed in the manifest and run by `run_pipeline.R` but missing from v0.1.0. It writes `methods.md` from the run's recorded parameters, summaries and package versions, then renders `templates/Report.Rmd` to `report.html`. A missing pandoc or failed render leaves `methods.md` in place and logs a warning.
- The report header now shows the code-bank commit (passed in by step 07, since `provenance.json` is written after the last step).
- New example dataset `gse60450` (mouse mammary gland, luminal vs basal; GEO GSE60450) in `make_example_data.R`, with `config/examples/gse60450.yaml`.
- Step 01 ignores numeric gene annotation columns written by counting tools (`Length`, `Start`, `End`, ...) instead of failing on them as unknown samples (new check id `annotation_columns`).
- NCBI Entrez Gene IDs are accepted (e.g. GEO's NCBI-generated raw count tables). Human or mouse is chosen by looking the IDs up in org.Hs.eg.db / org.Mm.eg.db.

## v0.1.0

First version of the code bank.

- Steps 01-07: input validation, QC, exploratory analysis (VST, PCA, outlier flags), DESeq2 differential expression with apeglm shrinkage, DE figures, pathway enrichment (fgsea ORA + GSEA on MSigDB or custom GMT gene sets), HTML report with an auto-generated methods paragraph.
- `run_pipeline.R` with run folders, per-step logs, provenance and run summaries; `--based_on` for reruns with changed settings.
- `manifest/tools.yaml` describing every step, parameter and output for the agent.
- Tests on airway and pasilla, including a comparison against a hand-written DESeq2 analysis and planted input errors.
- Developed and tested with R 4.3.3, DESeq2 1.42.0, apeglm 1.27.0, fgsea 1.39.4, msigdbr 26.1.1.
