#!/usr/bin/env Rscript
# =============================================================================
# tests/check_manifest.R  --  is manifest/tools.yaml consistent? (run before every release)
# =============================================================================
# For every module under `modules:`:
#   - every key of its defaults file is a manifest parameter, and every parameter has a default
#   - its entry point, defaults file, tests and step scripts exist
#   - its `site` block uses only the building blocks the website knows (docs/ADDING_A_MODULE.md)
#     and only names parameters the module has
# The website runs the same checks (and more) when a release arrives.
# Exit code 0 = no problems.
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "tests/x"), ".."))
})
m <- yaml::read_yaml(file.path(REPO_ROOT, "manifest", "tools.yaml"))
problems <- character(0)
bad <- function(...) problems <<- c(problems, paste0(...))

HANDLERS <- c("file", "files", "sample_folders", "annotation")
FORM_TYPES <- c("number", "integer", "enum", "bool", "text", "run_group")
RUN_VALUES <- paste0("run.", c("organism", "dir", "condition", "reference_level", "test_level", "covariates", "metadata_file", "samples"))
COLUMN_FORMATS <- c("num0", "num2", "p", "flag"); COLUMN_STYLES <- c("small", "sign"); TILE_STYLES <- c("up", "down")
exists_in_repo <- function(p) file.exists(file.path(REPO_ROOT, p))
as_chr <- function(x) if (is.null(x)) character(0) else as.character(unlist(x))

for (name in names(m$modules)) {
  mod <- m$modules[[name]]
  site <- mod$site
  if (is.null(site)) { bad(name, ": no site block (attach_to: run, page or none)"); next }
  if (!identical(site$attach_to, "none")) {
    for (k in c("entry_point", "defaults_file", "tests")) {
      if (is.null(mod[[k]])) bad(name, ": no ", k) else if (!exists_in_repo(mod[[k]])) bad(name, ": ", k, " not found: ", mod[[k]])
    }
    sdir <- dirname(mod$entry_point %||% "x")
    for (s in mod$steps) if (!exists_in_repo(file.path(sdir, paste0(s$name, ".R")))) bad(name, ": step script not found: ", s$name, ".R")
  }
  params <- names(mod$parameters)
  if (!is.null(mod$defaults_file) && exists_in_repo(mod$defaults_file)) {
    d <- names(yaml::read_yaml(file.path(REPO_ROOT, mod$defaults_file)))
    for (k in setdiff(d, params)) bad(name, ": default '", k, "' is not listed under parameters")
    for (k in setdiff(params, d)) bad(name, ": parameter '", k, "' has no default in ", mod$defaults_file)
  }
  if (!site$attach_to %in% c("run", "page", "none")) bad(name, ": site.attach_to must be run, page or none")
  if (!identical(site$attach_to, "run")) next
  for (i in site$inputs) {
    if (!i$handler %in% HANDLERS) bad(name, ": input '", i$name, "' uses handler '", i$handler, "' (known: ", paste(HANDLERS, collapse = ", "), ")")
    for (p in as_chr(i$sets)) if (!p %in% params) bad(name, ": input '", i$name, "' sets '", p, "', not a parameter")
  }
  for (f in site$form) {
    if (!f$type %in% FORM_TYPES) bad(name, ": form field '", f$name, "' has type '", f$type, "'")
    for (p in c(f$name, f$pair)) if (!p %in% params) bad(name, ": form field '", p, "' is not a parameter")
  }
  for (k in names(site$from_run)) {
    if (!k %in% params) bad(name, ": from_run sets '", k, "', not a parameter")
    v <- site$from_run[[k]]
    if (is.character(v) && startsWith(v, "run.") && !v %in% RUN_VALUES) bad(name, ": from_run uses unknown '", v, "'")
  }
  for (k in c(names(site$fixed), as_chr(site$inherit))) if (!k %in% params) bad(name, ": '", k, "' (fixed / inherit) is not a parameter")
  res <- site$results
  srcs <- names(res$sources)
  for (t in res$tiles) {
    if (!sub("\\..*", "", t$value) %in% srcs) bad(name, ": tile '", t$label, "' reads '", t$value, "' but results.sources lacks it")
    if (!is.null(t$style) && !t$style %in% TILE_STYLES) bad(name, ": tile '", t$label, "' style '", t$style, "'")
  }
  for (t in res$top_tables) for (c in t$columns) {
    if (!is.null(c$format) && !identical(c$format, TRUE) && !c$format %in% COLUMN_FORMATS) bad(name, ": column '", c$label, "' format '", c$format, "'")
    if (isTRUE(c$format) || isFALSE(c$format)) bad(name, ": column '", c$label, "' format was read as true/false: write it in quotes")
    if (!is.null(c$style) && !c$style %in% COLUMN_STYLES) bad(name, ": column '", c$label, "' style '", c$style, "'")
  }
  for (k in names(site$rerun$data)) if (!k %in% params) bad(name, ": rerun.data names '", k, "', not a parameter")
}

if (length(problems)) {
  cat("Manifest problems:\n", paste0("  - ", problems, "\n"), sep = "")
  quit(save = "no", status = 1)
}
cat(sprintf("Manifest OK: %d modules checked.\n", length(m$modules)))
