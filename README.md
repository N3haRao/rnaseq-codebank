# RNA-seq Code Bank

A tested, reproducible collection of R scripts for bulk RNA-seq differential expression analysis, built to be run by an AI agent behind a no-code website for wet-lab biologists. It can also be used directly from the command line.

You give it a **raw counts table** and a **sample table**. It checks the files, runs quality control, PCA, DESeq2 differential expression, volcano/MA/heatmap figures and pathway enrichment, and writes an HTML report with a methods paragraph. Every run is saved in its own folder together with the exact code version, parameters and package versions used.

This repository is the first layer ("code bank") of the RNA-seq analysis agent described in the project proposal. The agent never writes analysis code: it only calls the scripts here, with parameters described in [`manifest/tools.yaml`](manifest/tools.yaml).

---

## Contents

1. [Quick start](#1-quick-start)
2. [Requirements and installation](#2-requirements-and-installation)
3. [Running the example](#3-running-the-example)
4. [Analysing your own data](#4-analysing-your-own-data)
5. [The analysis steps](#5-the-analysis-steps)
6. [Parameters](#6-parameters)
7. [Outputs](#7-outputs)
8. [Changing settings and rerunning](#8-changing-settings-and-rerunning)
9. [Reproducibility](#9-reproducibility)
10. [Testing](#10-testing)
11. [Repository layout](#11-repository-layout)
12. [How the agent uses this repository](#12-how-the-agent-uses-this-repository)
13. [Troubleshooting](#13-troubleshooting)
14. [Status, limitations and citations](#14-status-limitations-and-citations)

---

## 1. Quick start

From a terminal, in the repository folder:

```bash
Rscript scripts/00_setup_environment.R                     # install packages (once)
Rscript scripts/make_example_data.R                        # write the example datasets
Rscript scripts/run_pipeline.R --config config/examples/airway.yaml
```

The last command prints where the results went, for example `runs/20261019-111338_airway/`. Open `07_report/report.html` in that folder.

---

## 2. Requirements and installation

**You need:** R 4.3 or newer, internet access for the first installation, and about 3 GB of disk space for packages. [pandoc](https://pandoc.org/installing.html) is needed for the HTML report; it comes with RStudio, and on Linux it can be installed with `sudo apt install pandoc`. Everything else still works without it.

Choose one of the three ways below.

### Option A: install packages directly (simplest)

```bash
git clone https://github.com/N3haRao/rnaseq-codebank.git
cd rnaseq-codebank
Rscript scripts/00_setup_environment.R
```

The setup script installs whatever is missing from CRAN and Bioconductor, downloads the MSigDB gene-set database once (it is cached afterwards), and prints a table of every package and version. To only check what is installed: `Rscript scripts/00_setup_environment.R --check`.

On Linux, compiling some packages needs system libraries. On Ubuntu/Debian:

```bash
sudo apt install build-essential libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev libpng-dev libjpeg-dev pandoc
```

### Option B: renv (recommended for the lab: pins exact versions)

[renv](https://rstudio.github.io/renv/) gives the project its own package library and records exact versions in `renv.lock`, so everyone in the lab, and any future server, runs the same versions.

First time, by one person:

```bash
Rscript -e 'install.packages("renv"); renv::init(bare = TRUE)'
Rscript scripts/00_setup_environment.R
Rscript -e 'renv::snapshot(type = "all")'
git add renv.lock .Rprofile renv/activate.R renv/settings.json
git commit -m "Pin package versions with renv"
```

Everyone else, after cloning:

```bash
Rscript -e 'install.packages("renv"); renv::restore()'
```

### Option C: Docker

```bash
docker build -t rnaseq-codebank .
docker run --rm rnaseq-codebank Rscript tests/run_tests.R
```

See the comments in the [`Dockerfile`](Dockerfile) for running your own data in the container and for pinning a Bioconductor release.

### Packages used

| Purpose | Packages |
| --- | --- |
| Differential expression | DESeq2, apeglm (ashr optional) |
| Pathway enrichment | fgsea, msigdbr |
| Gene annotation | AnnotationDbi, org.Hs.eg.db, org.Mm.eg.db |
| Plots | ggplot2, ggrepel, pheatmap |
| Plumbing | optparse, yaml, jsonlite, data.table, readxl (optional) |
| Report | rmarkdown, knitr, pandoc |
| Example data | airway, pasilla |

---

## 3. Running the example

Three public datasets are used for examples and tests. Create their input files with:

```bash
Rscript scripts/make_example_data.R
```

| Dataset | Experiment | Design | What it tests |
| --- | --- | --- | --- |
| **airway** (Himes et al. 2014) | Human airway smooth muscle cells, dexamethasone vs control, 4 donors | `~ cell_line + condition` (paired) | Standard two-group analysis with a paired design; a well-known result |
| **pasilla** (Brooks et al. 2011) | Drosophila S2 cells, pasilla knock-down vs untreated | `~ type + condition` | A technical covariate; an organism without MSigDB gene sets |
| **gse60450** (Fu et al. 2015) | Mouse mammary gland, luminal vs basal cells from virgin, pregnant and lactating mice, 2 replicates each | `~ status + cell_type` | A file as downloaded from GEO: Entrez IDs and an extra `Length` column; mouse annotation and pathways |

Run airway:

```bash
Rscript scripts/run_pipeline.R --config config/examples/airway.yaml
```

The console shows each step and a summary:

```
[11:13:38] INFO  New run: .../runs/20261019-111338_airway
[11:13:38] INFO  Step 1/7: 01_validate_inputs ...
...
================================================================
 Run:     20261019-111338_airway
 Status:  SUCCESS
 Compare: dexamethasone vs control (column 'condition')
 Design:  ~ cell_line + condition
 DEGs:    760 significant (412 up, 348 down); padj < 0.05 and |shrunken log2FC| > 1
 Folder:  .../runs/20261019-111338_airway
 Report:  07_report/report.html
================================================================
```

The whole run takes about one minute on a laptop. The top genes (CACNB2, SPARCL1, DUSP1, PER1, KLF15, FKBP5, ZBTB16) are the known glucocorticoid-response genes for this dataset.

---

## 4. Analysing your own data

### Step 1: prepare two files

See [`docs/INPUT_FORMAT.md`](docs/INPUT_FORMAT.md) for details and examples. In short:

- **Counts table:** genes in rows, samples in columns, gene IDs in the first column, **raw counts** (not TPM/FPKM/log). Blank template: [`templates/counts_template.csv`](templates/counts_template.csv).
- **Sample table:** one row per sample, sample IDs in the first column (matching the counts column names exactly), and a column with the groups. Blank template: [`templates/metadata_template.csv`](templates/metadata_template.csv).

### Step 2: write a run config

Copy [`config/examples/template.yaml`](config/examples/template.yaml), for example to `config/my_experiment.yaml`, and fill it in:

```yaml
run_name: ko_vs_wt
counts_file: /path/to/counts.csv
metadata_file: /path/to/samples.csv
organism: mouse
condition: genotype          # column in the sample table holding the groups
reference_level: WT          # control group
test_level: KO               # compared against the control
covariates: [batch]          # optional; for paired samples put the donor/mouse column here
```

Only list what differs from the defaults; every option is explained in [`config/default_params.yaml`](config/default_params.yaml).

### Step 3: check the files first (optional but quick)

```bash
Rscript scripts/run_pipeline.R --config config/my_experiment.yaml --steps 1
```

Read `01_validate_inputs/validation_report.md` in the new run folder. Every problem is described in plain language, for example:

> FAIL | sample_matching | These sample columns in the counts file have no row in the metadata: 'KO_3' (closest metadata ID: 'KO3'). Sample names must match exactly between the two files.

### Step 4: run everything

```bash
Rscript scripts/run_pipeline.R --config config/my_experiment.yaml
```

### Without a config file

All the main settings can be given as flags:

```bash
Rscript scripts/run_pipeline.R \
  --counts_file counts.csv --metadata_file samples.csv \
  --condition genotype --reference_level WT --test_level KO \
  --covariates batch --organism mouse --run_name ko_vs_wt
```

Any other parameter can be set with `--set`, separating several with `;`:

```bash
Rscript scripts/run_pipeline.R --config config/my_experiment.yaml \
  --set "padj_cutoff=0.1;log2fc_cutoff=0.58;gene_set_collections=[H, REACTOME]"
```

---

## 5. The analysis steps

`run_pipeline.R` runs seven steps in order, each in its own R process with its own log. Steps 1-5 are **critical** (the pipeline stops if one fails). Steps 6-7 are not: if pathway analysis fails, for example because there is no internet connection to fetch gene sets, the differential expression results are still complete and valid.

| # | Script | What it does |
| --- | --- | --- |
| 1 | [`01_validate_inputs.R`](scripts/01_validate_inputs.R) | Checks 26 things about the files (raw counts vs TPM/log, missing values, sample names matching, replicates, confounded covariates, ...). Explains every problem. Fills in open settings (reference group, organism, filter threshold) and records them. |
| 2 | [`02_qc.R`](scripts/02_qc.R) | Library sizes, detected genes, expression distributions. Flags samples with low or very small libraries. |
| 3 | [`03_exploratory.R`](scripts/03_exploratory.R) | Variance-stabilising transformation, PCA (coloured by group and by each covariate), sample-distance heatmap. Flags samples that cluster with the wrong group (possible mislabel) or sit far from their own group. |
| 4 | [`04_differential_expression.R`](scripts/04_differential_expression.R) | DESeq2 with the chosen design; Wald test; Benjamini-Hochberg adjusted p-values; apeglm fold-change shrinkage; gene symbols and names. |
| 5 | [`05_de_plots.R`](scripts/05_de_plots.R) | Volcano plot, MA plot, heatmap of top genes, p-value histogram. |
| 6 | [`06_enrichment.R`](scripts/06_enrichment.R) | Over-representation analysis of up and down genes, and pre-ranked GSEA, on MSigDB Hallmark, GO:BP, Reactome and KEGG (or your own `.gmt` files). |
| 7 | [`07_report.R`](scripts/07_report.R) | Methods paragraph written from what the run recorded, and a single-file HTML report. |

The statistics and the reasons for each choice are explained in plain language in [`docs/METHODS.md`](docs/METHODS.md).

Each step can also be run on its own on an existing run folder, which is handy when developing:

```bash
Rscript scripts/04_differential_expression.R --run_dir runs/20261019-111338_airway --padj_cutoff 0.1
Rscript scripts/04_differential_expression.R --help      # every option of a step
```

Note that running a step directly like this changes files in that run folder. For analyses you want to keep, use `run_pipeline.R --based_on` instead (section 8).

---

## 6. Parameters

All defaults live in [`config/default_params.yaml`](config/default_params.yaml), with a comment on each. The ones you are most likely to change:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `condition` | `condition` | Sample-table column holding the groups |
| `reference_level` | guessed | Control group. If not given, a control-like name (control, WT, vehicle, untreated, ...) is chosen and flagged for you to confirm |
| `test_level` | the other group | Group compared with the reference. Required when there are more than two groups |
| `covariates` | `[]` | Variables to adjust for, e.g. `[batch]`, `[sex, batch]`, or the donor column for paired samples |
| `exclude_samples` | `[]` | Sample IDs to leave out, e.g. `[KO_3]` |
| `organism` | `auto` | `human`, `mouse`, or `other` (no gene names or MSigDB pathways). Detected from the gene IDs |
| `padj_cutoff` | `0.05` | Adjusted p-value threshold |
| `log2fc_cutoff` | `1` | Absolute log2 fold-change threshold (1 = 2-fold; 0.58 = 1.5-fold), on the shrunken fold change |
| `lfc_shrinkage` | `apeglm` | `apeglm`, `ashr`, `normal` or `none` |
| `min_count`, `min_samples` | `10`, smaller group size | Low-count filter: keep genes with at least `min_count` reads in at least `min_samples` samples |
| `gene_set_collections` | `[H, GO:BP, REACTOME, KEGG]` | MSigDB collections for pathway analysis (also `GO:MF`, `GO:CC`, `WIKIPATHWAYS`) |
| `gene_set_source` / `gmt_files` | `msigdbr` / `[]` | Use `gmt` or `both` to add your own gene sets |
| `run_enrichment` | `true` | Set `false` to skip pathway analysis |

Parameters are combined in this order, later ones winning: defaults → run config (`--config`) → an earlier run (`--based_on`) → flags and `--set`. Step 01 then fills in anything left open and writes the final values to `01_validate_inputs/resolved_params.yaml`, which every later step reads.

---

## 7. Outputs

Each run creates `runs/<date-time>_<run_name>/`. The most useful files:

| File | What it is |
| --- | --- |
| `07_report/report.html` | Everything in one page: summary, warnings, all figures, top genes and pathways, methods, software versions |
| `run_summary.json` | Headline numbers and every warning in one small file (what the chat agent reads first) |
| `04_differential_expression/de_results_all.csv` | Every tested gene: fold change, p-value, padj, significance, direction, mean counts per group |
| `04_differential_expression/de_results_significant.csv` | Significant genes only |
| `05_de_plots/volcano_plot.png` | Volcano plot |
| `03_exploratory/pca_plot.png` | PCA |
| `06_enrichment/<collection>/gsea.csv` | GSEA results per collection |
| `07_report/methods.md` | Methods paragraph for a paper |
| `provenance.json` | Exact code commit, package versions, parameters and input checksums |

Every file and table column is described in [`docs/OUTPUTS.md`](docs/OUTPUTS.md).

---

## 8. Changing settings and rerunning

A finished run is never modified. To try different settings, start a new run from the old one's parameters and change only what you need:

```bash
Rscript scripts/run_pipeline.R --based_on runs/20261019-111338_airway --set "log2fc_cutoff=0.58;run_name=airway_fc1.5"
```

To leave a sample out (for example one flagged as an outlier, after checking the lab records):

```bash
Rscript scripts/run_pipeline.R --based_on runs/20261019-111338_airway --set "exclude_samples=[SRR1039508];run_name=airway_no508"
```

To finish a run that stopped (for example, pathway analysis failed while offline), resume it with the same parameters:

```bash
Rscript scripts/run_pipeline.R --run_dir runs/20261019-111338_airway --steps 6-7
```

---

## 9. Reproducibility

Every run records:

- the **git commit** of this repository (and a warning if there were uncommitted changes),
- **R and package versions** used by each step,
- all **parameters**, both requested and resolved,
- **MD5 checksums** of the input files, which are also copied into the run folder,
- timings and logs of each step.

To make results reproducible over time:

1. Commit before analysing, and tag releases: `git tag v0.1.0 && git push --tags`.
2. Pin package versions with renv (section 2, option B) or a fixed Docker tag.
3. Keep `provenance.json` with any results you share.

The test suite checks that a rerun from saved parameters gives an identical results table.

---

## 10. Testing

```bash
Rscript tests/run_tests.R                    # full suite, about 3 minutes
Rscript tests/run_tests.R --skip_enrichment  # when MSigDB cannot be downloaded
```

The suite checks that:

| Area | Check |
| --- | --- |
| Correctness | airway DE results (p-values, padj, fold changes, gene count) are identical to a hand-written DESeq2 analysis with the same settings |
| Biology | DUSP1, KLF15, PER1, TSC22D3, FKBP5 and ZBTB16 come out significant and up with dexamethasone |
| Reproducibility | a rerun with `--based_on` gives an identical table |
| Designs | pasilla runs with its covariate; enrichment is skipped cleanly for Drosophila |
| Safeguards | step 01 stops, with the right explanation, for: mismatched sample names, TPM values, log values, a missing value, a negative value, a missing group column, a batch confounded with the group, one replicate in a group, a misspelled reference group |
| Mislabels | a swapped control/treated label pair is flagged by the PCA check |
| Manifest | every parameter in `manifest/tools.yaml` has a default and every step has a script |

Run the tests after any change to the scripts, and before tagging a release. Output goes to `tests/output/` (not committed).

---

## 11. Repository layout

```
rnaseq-codebank/
├── README.md
├── CHANGELOG.md
├── Dockerfile
├── R/
│   └── utils.R                  shared helpers: parameters, logging, I/O, annotation, provenance, plot style
├── scripts/
│   ├── 00_setup_environment.R   install / check packages
│   ├── 01_validate_inputs.R
│   ├── 02_qc.R
│   ├── 03_exploratory.R
│   ├── 04_differential_expression.R
│   ├── 05_de_plots.R
│   ├── 06_enrichment.R
│   ├── 07_report.R
│   ├── run_pipeline.R           runs the steps as one tracked run
│   └── make_example_data.R      writes the airway / pasilla example inputs
├── config/
│   ├── default_params.yaml      every parameter, with defaults and explanations
│   └── examples/                airway.yaml, pasilla.yaml, template.yaml
├── manifest/
│   └── tools.yaml               machine-readable description of steps, parameters and outputs (for the agent)
├── templates/
│   ├── report.Rmd               HTML report template
│   ├── counts_template.csv      blank input templates for users
│   └── metadata_template.csv
├── docs/
│   ├── INPUT_FORMAT.md          how to prepare input files (for biologists)
│   ├── METHODS.md               statistical methods in plain language
│   ├── OUTPUTS.md               every output file and column
│   └── ADDING_A_STEP.md         how to add a new analysis
├── tests/
│   └── run_tests.R
├── data/example/                (generated, not committed)
└── runs/                        (analysis results, not committed)
```

---

## 12. How the agent uses this repository

The website and agent (later phases of the project) sit on top of this repository and follow fixed rules:

1. The agent reads [`manifest/tools.yaml`](manifest/tools.yaml) to learn which steps exist, their parameters (type, default, allowed values) and their outputs.
2. When a user uploads files, the agent runs step 1 only (`--steps 1`), reads `validation_report.json`, and proposes the comparison and covariates to the user in plain language. It waits for confirmation.
3. It then runs the full pipeline with the confirmed parameters. It never edits or writes analysis code.
4. To answer questions, it reads the files listed under `answers` in the manifest: `run_summary.json` first, then the specific tables. Numbers in its answers come from these files, never from memory.
5. "Rerun with X" becomes a new run with `--based_on <run> --set "X"`, so the original stays intact.

All scripts communicate through files in the run folder and exit codes, so the agent can be written in R (ellmer) or Python without changing anything here.

---

## 13. Troubleshooting

| Problem | What to do |
| --- | --- |
| A run stops at step 1 | Open `01_validate_inputs/validation_report.md` in the run folder; it says what is wrong and how to fix it |
| Any step fails | Read the last lines printed on screen, or the full log in `logs/<step>.log` |
| `there is no package called ...` | Run `Rscript scripts/00_setup_environment.R` |
| Enrichment fails with a download error | msigdbr fetches MSigDB on first use. Run the setup script once with internet access, or use `--set "run_enrichment=false"` |
| msigdbr asks to install `msigdbdf` | Accept: it is the data package that holds the gene sets in newer msigdbr versions |
| No `report.html` | pandoc is missing. Install it (or run from RStudio). `methods.md` and all other results are still there |
| Gene symbols are empty | The organism is `other`, the IDs were not recognised, or the annotation package (org.Hs.eg.db / org.Mm.eg.db) is not installed |
| "Fold-change shrinkage 'apeglm' failed" | apeglm is not installed; the run fell back to `normal` shrinkage and says so in the report. Install apeglm and rerun |
| Results differ between computers | Compare the package versions in the two `provenance.json` files; pin versions with renv |

---

## 14. Status, limitations and citations

**Version 0.1.0, trial stage, public data only.** Developed and tested with R 4.3.3, DESeq2 1.42.0, apeglm 1.27.0, fgsea 1.39.4 and msigdbr 26.1.1.

Current scope (v1): one two-group comparison per run, from a gene-level raw counts matrix, for human and mouse (other organisms get DE results without annotation or MSigDB pathways). Not included: FASTQ processing, time courses or interaction designs, single-cell or spatial data.

If you use these results in a publication, cite the underlying tools (listed in each run's `07_report/methods.md`):

- Love MI, Huber W, Anders S (2014). DESeq2. *Genome Biology* 15:550.
- Zhu A, Ibrahim JG, Love MI (2019). apeglm. *Bioinformatics* 35:2084–2092.
- Korotkevich G, et al. (2021). fgsea. *bioRxiv* 060012.
- Liberzon A, et al. (2015). MSigDB hallmark gene sets. *Cell Systems* 1:417–425.
