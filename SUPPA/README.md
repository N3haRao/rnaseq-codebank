# SUPPA2 module: alternative splicing from transcript abundance

SUPPA2 (Trincado et al. 2018) finds alternative-splicing events in the annotation and computes each event's
percent spliced in (PSI) from transcript abundance, so it works from Salmon output and needs no BAM files. It
complements the rMATS module (`splicing/`), which counts reads at the splice junctions; on the same samples the two
usually agree well on the strong events, and step SU04 shows how well.

| Event | Code | Meaning |
|---|---|---|
| Skipped exon | SE | an exon included or skipped |
| Alternative 5' / 3' splice site | A5 / A3 (from `SS`) | a longer or shorter exon at one end |
| Mutually exclusive exons | MX | one of two exons |
| Retained intron | RI | an intron kept in the mRNA |
| Alternative first / last exon | AF / AL (from `FL`) | a different first or last exon |

## Run it

```bash
Rscript suppa/scripts/run_suppa.R --quant_dir results/salmon --metadata_file samples.csv \
    --gtf_file ref/gencode.v44.primary_assembly.annotation.gtf --condition genotype \
    --reference_level WT --test_level KO
# with rMATS of the same comparison (rMATS --b1 = KO):
Rscript suppa/scripts/run_suppa.R ... --rmats_dir results/rmats/WT_vs_KO/tables
Rscript suppa/scripts/run_suppa.R --based_on runs/<run> --set "dpsi_cutoff=0.2"
```

**Input:** Salmon output (`<sample>/quant.sf`, from [`later/scripts/U02_salmon_quant.R`](../later/README.md)); the
sample table (first column = sample IDs = folder names); the GTF the Salmon index was built from (SUPPA finds the
events in it; GENCODE transcript names like `ENST...|ENSG...|...` are handled).

**SUPPA2 itself:** `suppa/install/install_suppa.sh` makes a conda environment at `$COPILOT_HOME/envs/suppa`, where
the module finds it (or set `suppa_command`). SUPPA2 is MIT-licensed. It needs older Python libraries (pandas < 2,
statsmodels 0.13); the script pins them.

## Steps

| Step | What it does | Main outputs |
|---|---|---|
| SU01 validate | SUPPA2 found and running; sample table, groups, replicates; quant.sf per sample; transcripts match the GTF (with or without versions); transcript TPM per group in SUPPA's format (every GTF transcript, 0 where Salmon has none) | `validation_report.json`, `<group>.tpm`, `resolved_params.yaml` |
| SU02 events + PSI | `suppa.py generateEvents -f ioe -e SE SS MX RI FL` (strict boundaries), merged; `psiPerEvent` per group | `events_all.ioe`, `<group>.psi`, `events_summary.json` |
| SU03 diffSplice | `suppa.py diffSplice -m empirical` with the reference first (dPSI = test − reference), `-gc` gene-level correction (default), `-pa` only for paired samples; significant events, genes, figures | `suppa_events_all.csv`, `suppa_events_significant.csv`, `suppa_genes_significant.csv`, `suppa_summary.json`, `event_counts.png`, `dpsi_volcano.png`, `psi_scatter.png` |
| SU04 vs rMATS | Optional (`rmats_dir`): events matched by event type, chromosome, strand and splice junctions; agreement of dPSI (Spearman) and of significant events, per type | `comparison_summary.json`, `comparison_by_type.csv`, `significant_in_both.csv`, `suppa_vs_rmats_dpsi.png` |
| SU05 report | Methods paragraph with references | `methods.md` |

**Significance:** p ≤ `pvalue_cutoff` (0.05; gene-corrected) and |dPSI| ≥ `dpsi_cutoff` (0.1). The lab's earlier
summaries used |dPSI| ≥ 0.2 for gene lists: `--set "dpsi_cutoff=0.2"`.

**Sign:** dPSI = PSI(test) − PSI(reference); positive = more inclusion in the test group, as rMATS with `--b1` = test
and MAJIQ with the reference as `-psi1`.

**Matching with rMATS (SU04)**, ported from the lab's `compare_SUPPA_rMATS.py`: a junction is "last base of the
upstream exon – first base of the downstream exon" (1-based; rMATS's `*Start_0base` + 1). SE/A5/A3/MX/RI ↔
SE/A5SS/A3SS/MXE/RI; RI on its retained intron only. For MX, SUPPA's inclusion form is the genomically first exon on
both strands while rMATS's is the 2nd exon on the − strand, so SUPPA's dPSI is flipped there (the Spearman for MX on
the − strand is reported separately to check this). `rmats_b1_is_test: false` if rMATS's `--b1` was the reference.

## Tests

```bash
Rscript suppa/tests/run_suppa_tests.R      # needs SUPPA2 (install_suppa.sh or SUPPA_COMMAND); skipped without it
```

19 tests on a synthetic annotation (300 cassette exons, 100 retained introns) and Salmon output (3 WT, 3 KO) with
planted changes: events found per type, planted events recovered with the right sign and few false positives,
summary and figures, an rMATS table matched by junctions, a `--based_on` rerun, and failures for a missing sample,
a non-matching GTF, missing SUPPA2 and an unknown group.
