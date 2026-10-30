# Changelog

Each release is tagged in git (`git tag v0.1.0`). Every run records the commit it used in `provenance.json`.

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
