# Statistical methods and the reasons behind them

This file explains each analysis choice in plain language. It is written for wet-lab users who want to understand their results, and as the reference the chat agent uses to answer "why" and "what does this mean" questions.

## Filtering low-count genes

Genes with fewer than `min_count` reads (default 10) in fewer than `min_samples` samples (default: the size of the smaller compared group) are removed before testing. A gene seen in only one or two samples at a handful of reads cannot show a reliable difference, and testing it only adds to the multiple-testing burden. Tying `min_samples` to the group size means a gene expressed in just one group is kept.

## Normalisation

DESeq2 estimates a **size factor** per sample with the median-of-ratios method. It corrects for sequencing depth and for differences in RNA composition, so that counts are comparable between samples. Normalised counts are the raw counts divided by the size factor (`04_differential_expression/normalized_counts.csv`). They are good for plotting a single gene, but the statistical test always works on the raw counts.

## Variance-stabilising transformation (VST)

Raw counts have a variance that grows with the mean, which makes highly expressed genes dominate distances and PCA. The VST puts counts on a log2-like scale where variance is roughly constant. It is used **only for visualisation** (PCA, sample distances, heatmaps), with `blind = TRUE` so the group labels do not influence the quality check.

## PCA and outlier flags

PCA summarises the 500 most variable genes into a few axes. Samples that are similar sit close together. When the groups separate along PC1 or PC2, the treatment effect is one of the largest sources of variation. `exploratory_summary.json` reports how much of PC1 and PC2 each metadata variable explains (R²), which helps spot batch effects.

Two automatic flags are raised, never acted on automatically:

- **clusters_with_other_group**: a sample sits closer to the centre of another group than to its own (the own-group centre is computed without the sample). Possible causes: a swapped label, a failed treatment, or a batch effect stronger than the treatment.
- **far_from_own_group** (6 or more samples only): a sample is unusually far from its own group (robust z > 3.5).

A flagged sample should be checked against lab records before anyone decides to remove it.

## The DESeq2 model

Counts for each gene are modelled with a **negative binomial** distribution, which handles the extra variability of biological replicates better than a Poisson or normal model. The design formula (e.g. `~ batch + condition`) lists what explains the counts; the last term is the comparison of interest and the others are adjusted for. DESeq2 shares information across genes to estimate each gene's variability (dispersion), which makes small experiments more stable. `dispersion_plot.png` shows this fit: black dots (per-gene) should scatter around the red line (trend).

## Testing and p-values

The **Wald test** asks whether the log2 fold change for the comparison differs from zero.

- **pvalue**: the chance of a difference at least this large if the gene were truly unchanged.
- **padj**: the p-value adjusted for testing thousands of genes at once (Benjamini-Hochberg). Keeping genes with padj < 0.05 means about 5% of the genes on the list are expected to be false positives. Use padj, not pvalue, to decide significance.
- **Independent filtering**: DESeq2 automatically leaves very low-count genes out of the multiple-testing correction, which increases power. Those genes get `padj = NA`; the `note` column says so.
- **Cook's distance**: a gene whose result is driven by one extreme sample gets `pvalue = NA` and a note.

## Fold changes and shrinkage

**log2FoldChange** compares the test group with the reference: +1 means twice as high in the test group, −1 means half, 0 means no change. Positive is always "higher in the test group".

Genes with low counts or high variability can show huge but unreliable fold changes. **Shrinkage** (apeglm by default) pulls those noisy estimates towards zero while barely changing well-measured genes. That makes the fold change more trustworthy for ranking, plots, and the fold-change cutoff. The unshrunken value is kept in `log2FoldChange_unshrunken`. P-values are not affected by shrinkage.

## What counts as "significant"

A gene is significant when **padj < 0.05 and |shrunken log2FC| > 1** (at least a 2-fold change). Both cutoffs can be changed; for subtle effects people often use |log2FC| > 0.58 (1.5-fold). A gene that is not significant is not proven unchanged; the experiment may simply lack the power to detect a small change.

## Pathway analysis

Two complementary approaches, both from the fgsea package:

**Over-representation analysis (ORA)** takes the list of significant up-regulated genes (and separately the down-regulated ones) and asks which gene sets contain more of them than expected by chance. The background is every gene that was tested, not the whole genome, which avoids inflated results. `fold_enrichment` is how many times more often than expected the set's genes appear in the list.

**Gene set enrichment analysis (GSEA)** uses every tested gene ranked by the DESeq2 test statistic, so no cutoff is needed. It detects pathways where many genes shift a little in the same direction, which ORA can miss. **NES** (normalised enrichment score) is positive when the set is shifted towards the test group and negative when shifted towards the reference. **Leading-edge genes** are the genes driving the result.

Gene sets come from MSigDB through the msigdbr package: Hallmark (50 well-defined biological states, a good first look), GO biological process, Reactome and KEGG. Mouse analyses use the human MSigDB with ortholog mapping. Gene sets with fewer than 15 or more than 500 tested genes are skipped, and p-values are adjusted per collection.

Enrichment results suggest biology; they do not prove that a pathway is active. Overlapping gene sets (common in GO) often appear together for the same few genes, so check the overlap and leading-edge genes.

## Reproducibility

Every run records the code-bank git commit, R and package versions, all parameters, input file checksums and timings in `provenance.json`. Rerunning from those parameters with the same code and package versions gives identical results (checked by `tests/run_tests.R`). Exact package versions are pinned with `renv.lock` (see the README).

## References

- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550.
- Zhu A, Ibrahim JG, Love MI (2019). Heavy-tailed prior distributions for sequence count data. *Bioinformatics* 35:2084–2092.
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate. *J R Stat Soc B* 57:289–300.
- Korotkevich G, et al. (2021). Fast gene set enrichment analysis. *bioRxiv* 060012.
- Subramanian A, et al. (2005). Gene set enrichment analysis. *PNAS* 102:15545–15550.
- Liberzon A, et al. (2015). The MSigDB hallmark gene set collection. *Cell Systems* 1:417–425.
