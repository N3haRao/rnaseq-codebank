# Splicing module: differential alternative splicing from rMATS output

This module starts from **rMATS output tables**, the same way the main pipeline starts from a counts table. It works on rMATS results made anywhere: on a cluster, by a core facility, or by [`later/scripts/U04_rmats.R`](../later/README.md). It uses the same conventions as the main pipeline: each run gets a folder with a parameter file, a log per step, provenance, a `run_summary.json`, and an entry under `modules` in [`manifest/tools.yaml`](../manifest/tools.yaml).

## Run it

```bash
Rscript splicing/scripts/run_splicing.R --rmats_dir path/to/rmats_output \
    --group_1_name KO --group_2_name WT --run_name ko_vs_wt_splicing

# with sample names for the plots, and other settings
Rscript splicing/scripts/run_splicing.R --rmats_dir path/to/rmats_output --group_1_name KO --group_2_name WT \
    --set "group_1_samples=[KO1, KO2, KO3];group_2_samples=[WT1, WT2, WT3];delta_psi_cutoff=0.2"

# from a config (U04 writes one per comparison), or a rerun with a change
Rscript splicing/scripts/run_splicing.R --config results/rmats/WT_vs_KO/splicing_config.yaml
Rscript splicing/scripts/run_splicing.R --based_on runs/<run> --set "fdr_cutoff=0.01"
```

**Input.** The folder rMATS wrote with `--od`, containing `SE.MATS.JCEC.txt`, `MXE.MATS.JCEC.txt`, `A3SS.MATS.JCEC.txt`, `A5SS.MATS.JCEC.txt` and `RI.MATS.JCEC.txt`. Use `count_type: JC` for the junction-only tables. rMATS must have been run with statistics (not `--statoff`).

**Groups.** Group 1 is rMATS `--b1` and group 2 is `--b2`. dPSI = PSI(group 1) − PSI(group 2), so a positive value ("up") means more inclusion of the alternative exon or region in group 1. U04 always puts the test group in b1 and the control in b2.

## Steps

| Step | What it does | Main outputs |
|---|---|---|
| S01 validate | Finds the tables, checks columns, that statistics are present, replicates and group names; standardises all events (1-based coordinates, mean PSI and reads per group) | `validation_report.json`, `events_validated.csv`, `resolved_params.yaml` |
| S02 results | Calls an event significant if FDR ≤ `fdr_cutoff`, \|dPSI\| ≥ `delta_psi_cutoff`, and there are at least `min_reads_per_sample` reads per sample in **both** groups | `splicing_events_all.csv`, `splicing_events_significant.csv`, `splicing_summary.json`, `event_counts.png`, `dpsi_volcano.png`, `psi_scatter.png`, `top_events_heatmap.png` |
| S03 enrichment | `fgsea::fora` on genes with significant events (all event types together, and per type), against all genes with an event that passed the read filter. Same gene sets as step 06 (MSigDB or GMT) | `<collection>/ora_*.csv`, dot plots, `enrichment_summary.json` |
| S04 report | Methods paragraph from the recorded parameters and results | `methods.md` |

Settings (all in [`config/default_params.yaml`](config/default_params.yaml)): `fdr_cutoff` 0.05, `delta_psi_cutoff` 0.1, `min_reads_per_sample` 10, `event_types`, `count_type`, `organism` (`auto` reads it from the gene IDs), plot colours and labels, and the gene-set settings of the main pipeline.

## Methods and choices

- **Statistics are rMATS's own.** The likelihood-ratio test on inclusion levels, with Benjamini-Hochberg FDR within each event type, is used as rMATS reports it. The read filter decides which events can be called significant; it does not recompute the FDR.
- **No GSEA for splicing.** A gene can have several events in opposite directions, so there is no single gene-level score to rank genes by. Only gene lists are tested (over-representation), using the genes that could have been detected as the background.
- **Changes from the lab's Python scripts.** The pathway step used to fill in random fold changes and p-values for the spliced genes; it now tests the gene lists directly. A run that fails stops with an error and a message instead of finishing silently.

## Tests

```bash
Rscript splicing/tests/run_splicing_tests.R
```

16 tests on synthetic rMATS tables with planted events. They check that exactly the planted events are found in the right direction, that low-coverage events are not called, that the planted pathway comes out of the enrichment, that reruns are reproducible, and that bad input (missing tables, `--statoff` output, missing columns, wrong count type) fails with the right message. They run offline.
