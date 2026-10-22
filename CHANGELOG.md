# Changelog

Each release is tagged in git (`git tag v0.1.0`). Every run records the commit it used in `provenance.json`.

## Unreleased

- Added `scripts/07_report.R`, listed in the manifest and run by `run_pipeline.R` but missing from v0.1.0. It writes `methods.md` from the run's recorded parameters, summaries and package versions, then renders `templates/Report.Rmd` to `report.html`. A missing pandoc or failed render leaves `methods.md` in place and logs a warning.
- The report header now shows the code-bank commit (passed in by step 07, since `provenance.json` is written after the last step).
- New example dataset `gse60450` (mouse mammary gland, luminal vs basal; GEO GSE60450) in `make_example_data.R`, with `config/examples/gse60450.yaml`.
- Step 01 ignores numeric gene annotation columns written by counting tools (`Length`, `Start`, `End`, ...) instead of failing on them as unknown samples (new check id `annotation_columns`).
- NCBI Entrez Gene IDs are accepted (e.g. GEO's NCBI-generated raw count tables). Human or mouse is chosen by looking the IDs up in org.Hs.eg.db / org.Mm.eg.db.

## v0.1.0 (2026-10-19)

First version of the code bank.

- Steps 01-07: input validation, QC, exploratory analysis (VST, PCA, outlier flags), DESeq2 differential expression with apeglm shrinkage, DE figures, pathway enrichment (fgsea ORA + GSEA on MSigDB or custom GMT gene sets), HTML report with an auto-generated methods paragraph.
- `run_pipeline.R` with run folders, per-step logs, provenance and run summaries; `--based_on` for reruns with changed settings.
- `manifest/tools.yaml` describing every step, parameter and output for the agent.
- Tests on airway and pasilla, including a comparison against a hand-written DESeq2 analysis and planted input errors.
- Developed and tested with R 4.3.3, DESeq2 1.42.0, apeglm 1.27.0, fgsea 1.39.4, msigdbr 26.1.1.
