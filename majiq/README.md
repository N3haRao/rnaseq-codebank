# MAJIQ module: local splicing variations from aligned reads

MAJIQ (Vaquero-Garcia et al. 2016, 2023) builds a splice graph from the annotation and the reads, including junctions
that are not annotated, and describes splicing as **local splicing variations** (LSVs): all the junctions that leave
(or enter) one exon. For each junction it gives the posterior probability that its inclusion changes between the
groups by at least a set amount. `voila modulize` then sorts the changes into event types. It complements rMATS
(`splicing/`, annotated events, read counts) and SUPPA2 (`suppa/`, transcript abundance).

**Licence:** MAJIQ is free for non-commercial use by individuals at academic or not-for-profit institutions, under
its academic licence; it is not open source and is not built into the site image. Each server installs it with
`majiq/install/install_majiq.sh` after registering at https://majiq.biociphers.org/app_download/. The module finds it
at `$COPILOT_HOME/envs/majiq/bin` and the licence file at `$COPILOT_HOME/majiq_license*.lic` (or set `majiq_bin_dir`,
`majiq_license_file`).

## Run it

```bash
Rscript majiq/scripts/run_majiq.R --bam_dir results/star --metadata_file samples.csv \
    --gff3_file ref/gencode.v44.primary_assembly.annotation.gff3.gz --condition genotype \
    --reference_level WT --test_level KO
Rscript majiq/scripts/run_majiq.R --based_on runs/<run> --set "probability_threshold=0.9;threads=16"
```

**Input:**
- sorted, indexed BAMs: `<sample>.bam`, `<sample>.sorted.bam` or `<sample>/Aligned.sortedByCoord.out.bam` (STAR in
  [`later/`](../later/README.md)), each with its `.bai` (made with samtools when missing and samtools is there);
- the sample table (first column = sample IDs);
- the annotation as **GFF3** of the same release as the GTF used for alignment (MAJIQ does not read GTF), e.g.
  `gencode.v44.primary_assembly.annotation.gff3.gz` for GENCODE 44.

## Steps (the lab's MAJIQ v3 workflow)

| Step | Command | Main outputs |
|---|---|---|
| MJ01 validate | MAJIQ, voila (and majiq-v3), licence; sample table, groups; GFF3 (not a GTF); BAMs + index, linked as `bams/<sample>.bam`; `experiments.tsv` (path, group) | `validation_report.json`, `experiments.tsv`, `resolved_params.yaml` |
| MJ02 build | `majiq build --license … <gff3> experiments.tsv build --strandness AUTO --nthreads N` | `build/splicegraph.zarr`, `build/<sample>.sj` |
| MJ03 quantify | `majiq psi-coverage` per group; `majiq-v3 sg-coverage` per group (needed by modulize) | `psicov/<group>.psicov`, `sgc/<group>.sgc` |
| MJ04 deltapsi | `majiq deltapsi -psi1 <reference> -psi2 <test> -n <ref> <test> --changing-threshold 0.2`; `voila tsv` (all LSVs with read counts; changing LSVs at 0.2 / 0.95); `voila modulize … sgc/ --decomplexify-psi-threshold 0.05` | `voila_all.tsv`, `voila_changing.tsv`, `modulize/` |
| MJ05 results | Per LSV the junction with the highest P(\|dPSI\| ≥ 0.2); changing when ≥ 0.95; genes; modulize event types; figures | `majiq_lsvs_all.csv`, `majiq_lsvs_changing.csv`, `majiq_genes_changing.csv`, `modulize_event_types.csv`, `majiq_summary.json`, `lsv_volcano.png`, `event_types.png` |
| MJ06 report | Methods paragraph | `methods.md` |

**Sign:** dPSI = PSI(test) − PSI(reference) (reference as `-psi1`), the same as rMATS (`--b1` = test) and SUPPA2.

**voila's column names** differ between versions. MJ05 finds them by name (`mean_dpsi_per_lsv_junction`,
`probability_changing`, `<group>_mean_psi`, `lsv_id`, ...); if a needed column is missing it stops with the list of
the columns it found, and the voila tables of MJ04 are still there.

## Tests

```bash
Rscript majiq/tests/run_majiq_tests.R
```

24 tests with stand-in `majiq` / `majiq-v3` / `voila` programs (MAJIQ needs a licence and real BAMs): the commands
match the lab's workflow (licence, strandness, reference as `-psi1`, thresholds, modulize with the sgc folder), the
voila tables are read correctly (JSON comment block, `;`-separated junction values), changing LSVs and directions,
modulize counts, figures, `run_modulize: false`, unrecognised voila columns, and failures for missing MAJIQ, licence,
BAM, BAM index and a GTF given instead of a GFF3. With MAJIQ installed, `MAJIQ_TEST_REAL=1` also checks that it starts.
