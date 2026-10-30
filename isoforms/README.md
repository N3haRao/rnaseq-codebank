# Isoform module: differential isoform expression (DEI) and usage (DUI)

Gene-level analysis adds up all of a gene's transcripts. This module looks at the transcripts themselves. It asks two different questions:

| | Question | Method | Example |
|---|---|---|---|
| **DEI**, differential isoform expression | Is this transcript more or less abundant? | DESeq2 per transcript: the main pipeline's method, with the same filter, design, test and fold-change shrinkage | A transcript is 4× higher in KO |
| **DUI**, differential isoform usage | Does this transcript make up a larger or smaller share of its gene? (an isoform switch) | DEXSeq on transcript shares within each gene, with stageR two-stage FDR control (Love, Soneson & Patro 2018) | The gene's total is unchanged, but transcript A goes from 80% to 20% of it and transcript B takes over |

A gene that is simply more expressed shows up in DEI but **not** in DUI. A pure isoform switch shows up in DUI, and its transcripts usually in DEI too.

## Run it

```bash
Rscript isoforms/scripts/run_isoforms.R --quant_dir results/salmon --metadata_file samples.csv \
    --gtf_file ref/genes.gtf --condition genotype --reference_level WT --test_level KO \
    --set "covariates=[batch];threads=4"
Rscript isoforms/scripts/run_isoforms.R --based_on runs/<run> --set "min_delta_usage=0.2"
```

**Input:**
- **Quantification:** Salmon output with one folder per sample (`<sample>/quant.sf`, made by [`later/scripts/U02_salmon_quant.R`](../later/README.md)), or kallisto output (`quant_type: kallisto`).
- **Sample table:** its first column holds the sample IDs, which must match the folder names.
- **Annotation:** the GTF the transcriptome came from, or a `tx2gene_file`.

Salmon with bootstraps and bias correction (as U02 runs it) is recommended for DUI.

## Steps

| Step | What it does | Main outputs |
|---|---|---|
| I01 validate | Checks the sample table, groups and replicates, that the design can be estimated, the quantification files, and the transcript-to-gene match; imports with tximport (transcript level) | `validation_report.json`, `transcripts.csv`, `resolved_params.yaml` |
| I02 DEI | DESeq2 on transcript counts with effective-length offsets; low-count filter; Wald test; apeglm shrinkage | `dei_results_all.csv` (transcript, gene, log2FC, padj, mean TPM per group), `dei_results_significant.csv`, `dei_volcano.png`, `dei_ma_plot.png` |
| I03 DUI | Counts scaled for transcript usage (`dtuScaledTPM`); DRIMSeq-style filter; DEXSeq on the two compared groups; stageR (genes screened, then transcripts confirmed; overall FDR `dui_alpha`) | `dui_transcripts_all.csv` (share in each group, Δusage, gene and transcript adjusted p), `dui_genes.csv` (isoform switches), `dui_volcano.png`, `dui_top_genes_usage.png` |
| I04 report | Methods paragraph with references | `methods.md` |

**Significance rules:**
- **DEI:** padj < `padj_cutoff` and |shrunken log2FC| > `log2fc_cutoff`. These are the main pipeline's rules, with the same defaults (0.05 and 1).
- **DUI:** a transcript is called when its gene passes the stageR screen, the transcript is confirmed, and |Δusage| ≥ `min_delta_usage` (0.1, i.e. 10 percentage points of the gene).
- **Isoform switch:** a gene where one transcript gains share and another loses it.

## Tests

```bash
Rscript isoforms/tests/run_isoform_tests.R     # needs DESeq2, DEXSeq, stageR, tximport
```

15 tests on synthetic Salmon output for 800 genes. 25 genes are 4× up in KO with unchanged shares; 25 genes switch their main transcript from 80% to 20%. The tests check that:

- DEI finds the expression genes;
- DUI finds the switch genes and **not** the expression-only genes;
- the switching transcript loses about 0.6 of its share;
- false calls are rare;
- reruns are identical;
- bad input fails with the right message (missing quantification, a covariate confounded with the groups, annotation that doesn't match).
