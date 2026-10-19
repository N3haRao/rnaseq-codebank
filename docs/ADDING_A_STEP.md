# Adding a new analysis to the code bank

The agent can only run what is in the code bank, so new analyses are added here as tested scripts. Every script follows the same contract, which keeps the agent simple: it never needs to know how a script works, only what the manifest says about it.

## The contract

1. **One script, one step**, in `scripts/`, named `NN_short_name.R`. Its output folder in a run has the same name.
2. **Starts with the standard header**: a comment block saying what it does, its outputs and usage, then the three lines that find `REPO_ROOT` and source `R/utils.R` (copy them from any existing step).
3. **Takes `--run_dir`** plus `--config`, `--set` and one flag per parameter it reads. Parameters are loaded with `load_params(opt)`, which merges defaults, the run config, step 01's resolved values and the flags.
4. **Reads only from the run folder** (validated inputs, earlier steps' outputs). No hard-coded paths, no files outside the run.
5. **Writes only into `step_dir(run_dir, STEP)`**, and finishes with `write_step_provenance()` listing the parameters and packages it used.
6. **Writes a small JSON summary** with the headline numbers, so the agent can answer questions without reading big tables.
7. **Fails with `fail("plain-language message")`**, which exits with code 1.
8. **Is deterministic**: set a seed for anything random (`set.seed(params$seed)`).

## Checklist

- [ ] Script written following the contract above
- [ ] New parameters added to `config/default_params.yaml`, with a comment
- [ ] Entry added to `manifest/tools.yaml` (parameters with type and default, outputs, description)
- [ ] If it belongs in the standard pipeline: added to `PIPELINE_STEPS` in `R/utils.R` (and marked `critical` or not)
- [ ] If it should appear in the report: a section added to `templates/report.Rmd`, and a sentence to the methods text in `scripts/07_report.R`
- [ ] Tested on the airway example, with a check added to `tests/run_tests.R`
- [ ] `Rscript tests/run_tests.R` passes
- [ ] `CHANGELOG.md` updated and the version bumped in `R/utils.R` and `manifest/tools.yaml`

## Example: an alternative DE method

To add edgeR as an option, the cleanest route is a separate script `scripts/04b_edger.R` that writes the same table columns as step 04 into `04b_edger/`, plus a parameter `de_method: edger` that `run_pipeline.R` uses to pick which DE step to run. Steps 05 and 06 then read from whichever DE folder `resolved_params.yaml` names.
