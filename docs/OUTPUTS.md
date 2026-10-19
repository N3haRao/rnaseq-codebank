# What a run folder contains

Every run lives in its own folder, `runs/<YYYYmmdd-HHMMSS>_<run_name>/`. Nothing in an existing run is ever overwritten; a rerun with different settings creates a new folder.

```
runs/20261019-111338_airway/
├── run_summary.json              headline results and all warnings (start here)
├── provenance.json               code commit, R/package versions, parameters, input checksums, timings
├── params.yaml                   parameters requested for this run
├── inputs/                       copies of the uploaded counts and metadata
├── logs/                         full console output of each step
├── provenance/                   one record per step
├── 01_validate_inputs/
│   ├── validation_report.md      every input check, readable
│   ├── validation_report.json    same, for the agent
│   ├── resolved_params.yaml      final parameters (reference group, organism, filter, design)
│   ├── counts_validated.csv      cleaned integer counts
│   └── metadata_validated.csv    aligned sample table
├── 02_qc/
│   ├── qc_metrics.csv            per sample: reads, detected genes, % zeros, flags
│   ├── qc_summary.json
│   └── library_sizes.png, detected_genes.png, logcpm_boxplot.png, logcpm_density.png
├── 03_exploratory/
│   ├── pca_plot.png, pca_by_<covariate>.png, pca_scree.png, sample_distance_heatmap.png
│   ├── pca_coordinates.csv, pca_variance.csv
│   ├── outlier_flags.csv         possible mislabels / outliers
│   ├── exploratory_summary.json
│   └── vst_counts.csv.gz, size_factors.csv, vsd.rds
├── 04_differential_expression/
│   ├── de_results_all.csv        every tested gene (see columns below)
│   ├── de_results_significant.csv
│   ├── normalized_counts.csv
│   ├── de_summary.json           design, thresholds, counts, top genes
│   ├── dispersion_plot.png
│   └── dds.rds                   fitted DESeq2 object, for custom follow-up in R
├── 05_de_plots/
│   └── volcano_plot.png, ma_plot.png, heatmap_top_degs.png, pvalue_histogram.png
├── 06_enrichment/
│   ├── enrichment_summary.json
│   └── <collection>/             H, GO.BP, REACTOME, KEGG, ...
│       ├── ora_up.csv, ora_down.csv, gsea.csv
│       └── ora_up_dotplot.png, ora_down_dotplot.png, gsea_barplot.png
└── 07_report/
    ├── report.html               everything in one file (open in a browser; print to PDF if needed)
    └── methods.md                methods paragraph for a paper
```

## Columns of `de_results_all.csv`

| Column | Meaning |
| --- | --- |
| gene_id | ID from your counts file (Ensembl version removed) |
| gene_symbol, gene_name, entrez_id | Annotation (human and mouse only) |
| baseMean | Average normalised count across all samples |
| log2FoldChange | Shrunken log2 fold change, test vs reference (positive = higher in test) |
| lfcSE | Standard error of the log2 fold change |
| log2FoldChange_unshrunken | Fold change before shrinkage |
| stat | Wald statistic (used to rank genes for GSEA) |
| pvalue | Raw p-value |
| padj | Benjamini-Hochberg adjusted p-value |
| mean_normalized_<group> | Mean normalised count in each compared group |
| significant | TRUE when padj < cutoff and abs(log2FoldChange) > cutoff |
| direction | up, down, or not significant |
| note | Why pvalue or padj is missing, when it is |

## Pathway tables

`ora_up.csv` / `ora_down.csv`: pathway, pvalue, padj, n_overlap (significant genes in the set), set_size (tested genes in the set), fold_enrichment, overlap_genes.

`gsea.csv`: pathway, pvalue, padj, ES, NES (positive = higher in the test group), set_size, direction, leading_edge_genes.
