# Adding a module (and getting it onto the website without website code)

A module is a folder like `splicing/`, `isoforms/` or `compare/`: its own scripts, defaults, tests and run
folders, following the same rules as the main pipeline (parameters only; methods fixed by `methods:` in the
manifest). The website (the RNA-seq copilot) reads `manifest/tools.yaml`, so a module described there is
offered on the site **without changing the site's code**, as long as it only uses the building blocks below.

## 1. The module itself

```
mymodule/
  scripts/run_mymodule.R     entry point: --config <yaml> --runs_dir <folder>; makes runs/<time>_<run_name>/
  scripts/M01_....R ...      one script per step, each with its own log
  config/default_params.yaml every parameter with its default (run_name and seed included)
  R/mymodule_utils.R
  tests/run_mymodule_tests.R
  README.md
```

The entry point must:
- print `Step <i>/<n>: <step name> ...` when a step starts and `  done in <s> s` or `FAILED (exit code ...)` when
  it ends (the site shows progress from these lines);
- write `run_summary.json` with `status` = `success`, `success_with_errors` or `failed`;
- write a `validation_report.json` (in its first step's folder) with `checks: [{id, status, message}]`.

Copy `splicing/scripts/run_splicing.R`: it does all of this.

## 2. The manifest entry

Under `modules:` in `manifest/tools.yaml`: `description`, `entry_point`, `defaults_file`, `tests`, `steps`
(name, critical, description), `parameters` (**every** key of the defaults file), `outputs`, `answers`, and a
`site` block.

`Rscript tests/check_manifest.R` checks that parameters and defaults agree and that the `site` block only uses
known building blocks. Run it before every release.

## 3. The `site` block

```yaml
site:
  attach_to: run          # run  = offered on every finished main-pipeline run ("More analyses" tab)
                          # page = its own page on the site (needs site code; e.g. compare)
                          # none = not on the website (e.g. later/)
  title: My analysis
  description: >          # {test} {reference} {condition} {covariates} are filled in; **bold** works
    What it finds for {test} vs {reference}.
```

### Inputs (uploads)

Each input is stored and checked by a **handler**, which sets one module parameter:

| handler | what the user uploads | sets the parameter to | options |
|---|---|---|---|
| `file` | one file | its path | `accept` |
| `files` | files and/or .zip files | the folder holding the kept files | `keep` (regex on file names), `require` (regex at least one must match; `{form value}` allowed), `none_found`, `not_required_found` (messages) |
| `sample_folders` | one .zip with a folder per sample of the run | the folder (`<sample>/<file>`) | `file` (e.g. `quant.sf`); folder names must match the run's samples |
| `annotation` | optional GTF or tx2gene table | `gtf_file` or `tx2gene_file` | none; without an upload it builds tx2gene from GENCODE transcript names |

```yaml
  inputs:
    - {name: tables, handler: files, sets: my_dir, label: My tables, accept: [".txt", ".zip"], multiple: true,
       required: true, keep: '\.mytable\.txt$', help: Where the files come from.}
```

A new kind of upload needs a new handler on the site (a small Python function); everything else does not.

### Settings

```yaml
  from_run:                     # taken from the finished run
    organism: run.organism      # also: run.dir, run.condition, run.reference_level, run.test_level,
    de_run: run.dir             #       run.covariates, run.metadata_file (a copy of the validated sample table), run.samples
    run_name: "{test}_vs_{reference}_mine"
  fixed: {threads: 2}           # always this value
  inherit: [plot_palette, up_colour]   # copied from the run's advanced settings when set there
  form:                         # what the user can change; default and allowed values come from `parameters`
    - {name: cutoff, type: number, label: Cutoff, min: 0, max: 1, help: ...}
    - {name: method_variant, type: enum, labels: {a: "Variant A", b: "Variant B"}}
    - {name: draw_extra, type: bool}
    - {name: group_1_name, type: run_group, pair: group_2_name, default: test,
       notes: {reference: "Shown when the user picks the reference group first"}}
```

Form types: `number`, `integer`, `enum`, `bool`, `text`, `run_group` (the user picks which of the run's groups
is first; sets the field to that group's name and `pair` to the other one).

### Results

```yaml
  results:
    sources: {s: M02_results/summary.json}          # JSON files to read numbers from, by short name
    tiles: [{label: significant, value: s.n_significant, style: up}]      # style: up | down
    headline:                                        # sentences for the site, report and PowerPoint
      - "{s.n_significant} things were significant ({s.rule})."
      - {text: "Only shown when s.extra is set: {s.extra}.", when: s.extra}
      - {text: "Shown when s.status is not set.", unless: s.status}
    figures: [{path: M02_results/plot.png, title: My plot, caption: What it shows.}]   # shown if the file exists
    top_tables:
      - title: Top hits
        path: M02_results/significant.csv
        filter: {significant: "TRUE"}               # optional: keep rows with these values
        rows: 50
        empty: Nothing significant.
        columns:
          - {label: Gene, value: gene_name, fallback: gene_id}
          - {label: Where, template: "{chr}:{start}-{end}", style: small}
          - {label: log2 FC, value: log2FoldChange, format: num2, style: sign}   # format: num0 | num2 | p | flag
    tables: [M02_results/all.csv]                    # download buttons
    methods: M03_report/methods.md
    chat_files: [M02_results/summary.json, M02_results/all.csv]   # what the chat assistant is pointed to
  rerun:                                             # for the code .zip of a run
    data: {my_dir: tables}                           # parameters that point at input data -> folder name in the zip
    drop: [de_run]                                   # parameters not to keep
```

Write words like yes, no, on and off in quotes: YAML reads them as true and false.

## 4. Release

1. `Rscript mymodule/tests/run_mymodule_tests.R` and `Rscript tests/check_manifest.R` pass.
2. Bump `CODE_BANK_VERSION` in `R/utils.R` and `version:` in the manifest; add a CHANGELOG entry.
3. Commit and tag. The copilot's release agent picks it up: it checks the manifest against the site, runs the
   site's tests, makes any site changes needed in a staging copy, and asks you to approve.
