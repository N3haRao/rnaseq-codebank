# Compare module: do two studies agree?

This module compares the **results** of two to eight RNA-seq studies, for example your knockout against a published GEO dataset of the same knockout, or the same treatment in human and mouse cells.

It never merges raw counts. Counts from different studies differ in library preparation, sequencing depth, read length, annotation and batch, and a combined DESeq2 model would mostly measure those differences. So each study is analysed on its own, with its own normalisation and design, and the module compares what came out:

| Question | How it is answered | Step |
|---|---|---|
| Do the fold changes agree? | Pearson and Spearman correlation of log2 fold changes over the genes both studies tested, plus a scatter plot | C02 |
| Do the same genes come out significant, in the same direction? | Overlap of significant genes: up/up, down/down, up/down and down/up, each with a one-sided Fisher's exact test (the genes both studies tested are the universe), the Jaccard index and the share with the same direction | C02 |
| Do the rankings agree, without picking a cutoff? | Stratified rank-rank hypergeometric overlap (RRHO2): the overlap of the two ranked gene lists at every pair of thresholds, separately for same-direction and opposite-direction changes | C03 |
| Are the same pathways affected? | GSEA normalised enrichment scores (NES) across studies: a heatmap, the NES correlation, and the pathways significant in both studies | C04 |

## Run it

```bash
# two finished pipeline runs
Rscript compare/scripts/run_compare.R --studies "KO_ours=runs/20261101-120000_ko_vs_wt,GSE60450=runs/20261101-130000_gse60450"

# a config file, e.g. with a published table that compared the groups the other way round
Rscript compare/scripts/run_compare.R --config my_comparison.yaml

# rerun with other settings (the run folder keeps copies of every study's results)
Rscript compare/scripts/run_compare.R --based_on runs/<compare_run> --set "padj_cutoff=0.1;log2fc_cutoff=0.58"
```

`my_comparison.yaml`:

```yaml
run_name: ko_vs_published
studies:
  - {label: KO_ours, run_dir: runs/20261101-120000_ko_vs_wt}
  - {label: Smith2024, de_file: data/smith2024_supplementary_table3.xlsx, organism: mouse,
     comparison: "WT vs KO", flip: true}
```

**Studies:** each one is either a main-pipeline run folder or a DE results table (`de_file`).
- The table needs a gene column and a fold change plus an adjusted p value. DESeq2 (`log2FoldChange`, `padj`), edgeR (`logFC`, `FDR`) and limma (`logFC`, `adj.P.Val`) column names are recognised; `pvalue` and `stat` are used when present.
- For a table, the organism is detected from the gene IDs; set it when they are Entrez IDs.
- `flip: true` reverses a study's fold changes, for a study that compared WT vs KO when the others compared KO vs WT. If two studies' fold changes are negatively correlated, the run warns that one of them probably needs it.

## Choices the module makes, and how to change them

- **Matching genes** (`match_genes_by: auto`): same organism → the gene ID when every study uses the same type, else the Entrez ID, else the symbol. Human with mouse → mouse genes are mapped to their **one-to-one** human orthologs (babelgene, with at least `min_ortholog_support` of 12 databases agreeing); genes with no one-to-one ortholog are left out, and the number is reported.
- **One significance rule for every study** (`padj_cutoff` 0.05, `log2fc_cutoff` 1). Different cutoffs would make the overlaps meaningless. `use_study_cutoffs: true` keeps each run's own calls instead.
- **Ranking for RRHO** (`rank_by: auto`): the DESeq2 Wald statistic when every study has it, else sign(log2FC) × −log10(p) for all studies.
- **RRHO**: about 100 thresholds along each list (`rrho_step`), Benjamini-Yekutieli adjustment over the whole map (`rrho_correction`).
- **Pathways**: each run's own GSEA results are reused (step 06, fgsea on the Wald statistic). For a study given as a table, GSEA is run here the same way. Pathways are matched by name within each collection; MSigDB uses the same names for human and mouse. A study with no GSEA results (a table of an organism without gene sets) is left out of the pathway comparison; the others are still compared.

## Outputs

| Folder | Main files |
|---|---|
| `inputs/<study>/` | Copies of each study's results (DE table, settings, GSEA tables) so the comparison can be rerun on its own |
| `C01_harmonise_studies/` | `studies.json` (what each study compared, organism, genes matched, significant genes), `matched_genes.csv` (one row per gene, every study's log2FC, padj and direction), `validation_report.json` |
| `C02_compare_genes/` | Per pair, in `<a>_vs_<b>/`: `fc_scatter.png`, `overlap_matrix.png`, `genes_compared.csv` (each gene's category: both same direction / opposite / one study only / neither), `overlap_table.csv`. Across studies: `correlation_heatmap.png`, `genes_significant_in_several_studies.csv`, `intersections.png` (3+ studies), `gene_comparison.json` |
| `C03_rrho/` | Per pair: `rrho_map.png`, `rrho_grid.csv`, `rrho_best_genes_<quadrant>.csv` (the genes at the strongest cell); `rrho_summary.json` |
| `C04_compare_pathways/` | Per collection: `pathways_all_studies.csv`, `nes_heatmap.png`, and per pair `nes_scatter.png` and `pathways_compared.csv`; `pathway_comparison.json` |
| `C05_compare_report/` | `methods.md` |

## Reading the results

- **Correlation** of fold changes over all genes is often modest (0.2–0.5) even for studies that agree well, because most genes do not change and their fold changes are noise. Look at the correlation among genes significant in either study, the overlaps and RRHO too.
- **RRHO map**: a hot spot in the "up in both" corner (bottom left) and the "down in both" corner (top right) means the studies agree; hot spots in the other two corners mean they change genes in opposite directions.
- **Different organisms**: fewer genes are compared (only one-to-one orthologs), and species differences add to study differences.

## Tests

```bash
Rscript compare/tests/run_compare_tests.R     # needs DESeq2, fgsea, babelgene (about 2 minutes)
```

32 tests. Two synthetic studies are analysed with the real pipeline, and a third is an edgeR-style table that compared the groups the other way round. Planted genes are up in all studies, down in all, up in one study only, or up in one and down in the other. The tests check that:

- the planted genes land in the right overlap categories;
- the Fisher tests, correlations, RRHO quadrants and pathway comparison find the planted pattern;
- the flipped table agrees only after `flip: true`, with a warning without it;
- reruns are identical, and a stricter cutoff gives fewer genes;
- human and mouse tables are matched through orthologs;
- bad inputs fail with a clear message.

## References

- Plaisier SB, Taschereau R, Wong JA, Graeber TG. Rank-rank hypergeometric overlap: identification of statistically significant overlap between gene-expression signatures. *Nucleic Acids Res.* 2010;38:e169.
- Cahill KM, Huo Z, Tseng GC, Logan RW, Seney ML. Improved identification of concordant and discordant gene expression signatures using an updated rank-rank hypergeometric overlap approach. *Sci Rep.* 2018;8:9588.
