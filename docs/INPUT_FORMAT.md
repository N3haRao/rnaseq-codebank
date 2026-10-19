# Preparing your input files

You need two files: a **counts table** and a **sample table** (metadata). Blank templates are in `templates/`. Accepted formats are `.csv`, `.tsv`, `.txt` (tab or comma separated) and `.xlsx` (first sheet only).

## 1. Counts table

One row per gene, one column per sample. The first column holds the gene IDs.

| gene_id | WT_1 | WT_2 | WT_3 | KO_1 | KO_2 | KO_3 |
| --- | --- | --- | --- | --- | --- | --- |
| ENSMUSG00000000001 | 1523 | 1388 | 1610 | 1702 | 1590 | 1655 |
| ENSMUSG00000000003 | 0 | 0 | 0 | 0 | 0 | 0 |
| ENSMUSG00000000028 | 87 | 102 | 95 | 15 | 22 | 19 |

Requirements:

- **Raw read counts**, not TPM, FPKM, CPM, or log values. Raw counts are whole numbers and each sample adds up to millions of reads. They come from featureCounts, HTSeq-count, STAR `--quantMode GeneCounts`, or the `NumReads` / `expected_count` column of Salmon or RSEM (decimals there are fine; they are rounded).
- Gene IDs can be Ensembl IDs (`ENSG...`, `ENSMUSG...`, version suffixes like `.12` are removed) or gene symbols (`TP53`, `Trp53`).
- Use 0 for genes with no reads. Empty cells are not allowed.
- Extra text columns such as `gene_name` or `biotype` are ignored automatically.

If the gene IDs are not in the first column, set `gene_id_column` in your config.

## 2. Sample table (metadata)

One row per sample. The first column holds sample IDs that **exactly match** the column names of the counts table (same spelling, same upper/lower case).

| sample_id | genotype | sex | batch |
| --- | --- | --- | --- |
| WT_1 | WT | F | 1 |
| WT_2 | WT | M | 1 |
| WT_3 | WT | F | 2 |
| KO_1 | KO | M | 1 |
| KO_2 | KO | F | 2 |
| KO_3 | KO | M | 2 |

- One column holds the groups you want to compare (here `genotype`). Tell the pipeline its name with `condition`.
- Other columns are optional. Add any you may want to adjust for: batch, sex, sequencing run, donor.
- Every sample needs a value in the group column and in each covariate you use.
- The order of rows does not matter.

## Choosing the design

| Your experiment | Setting |
| --- | --- |
| Two groups, nothing else | `condition: genotype`, `covariates: []` |
| Samples processed in batches | `covariates: [batch]` |
| Paired samples (same donor / mouse / cell line in both groups) | `covariates: [donor]` |
| Mixed sexes | `covariates: [sex]` |
| More than two groups | Set `reference_level` and `test_level`; run once per comparison |

A covariate **cannot** be adjusted for if each of its values contains only one group (for example, all KO samples in batch 2 and all WT in batch 1). That is called confounding; the pipeline stops and explains it. No analysis can separate those two effects afterwards, so it is worth balancing groups across batches when planning experiments.

## Replicates

At least 2 biological replicates per compared group are required, and 3 or more are strongly recommended. With fewer, only large changes are detectable and a single unusual sample can drive the results.

## What happens if something is wrong

Step 01 checks the files before anything runs and writes `01_validate_inputs/validation_report.md`, which explains every problem in plain language (for example: "These sample columns in the counts file have no row in the metadata: 'KO_3' (closest metadata ID: 'KO3')").
