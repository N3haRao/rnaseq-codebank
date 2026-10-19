#!/usr/bin/env Rscript
# =============================================================================
# 00_setup_environment.R  --  install and check every package the code bank uses
# =============================================================================
# Run once from the repository root:
#   Rscript scripts/00_setup_environment.R            # install what is missing, then check
#   Rscript scripts/00_setup_environment.R --check    # only check, install nothing
#
# What it does
#   1. Installs BiocManager (if missing), then any missing CRAN and
#      Bioconductor packages listed below.
#   2. Downloads the MSigDB gene-set database used by msigdbr (first use only;
#      needs internet once, then it is cached on this computer).
#   3. Prints a table of every package and version, and whether pandoc (needed
#      for the HTML report) is available.
#
# Using renv (recommended, records exact versions in renv.lock):
#   Rscript -e 'install.packages("renv"); renv::init(bare = TRUE)'
#   Rscript scripts/00_setup_environment.R
#   Rscript -e 'renv::snapshot(type = "all")'
#   git add renv.lock .Rprofile renv/activate.R renv/settings.json
# Later, on another computer:  Rscript -e 'renv::restore()'
# =============================================================================

check_only <- "--check" %in% commandArgs(trailingOnly = TRUE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "optparse",   # command-line options for every script
  "yaml",       # parameter files
  "jsonlite",   # summaries and provenance records
  "data.table", # fast table reading/writing
  "ggplot2",    # plots
  "ggrepel",    # non-overlapping labels on plots
  "pheatmap",   # heatmaps
  "rmarkdown",  # HTML report
  "knitr",      # HTML report
  "readxl",     # optional: accept .xlsx uploads
  "msigdbr",    # MSigDB gene sets
  "ashr"        # optional: alternative fold-change shrinkage
)
bioc_pkgs <- c(
  "DESeq2",        # differential expression
  "apeglm",        # default fold-change shrinkage
  "fgsea",         # ORA and GSEA
  "AnnotationDbi", # gene annotation
  "org.Hs.eg.db",  # human gene annotation
  "org.Mm.eg.db",  # mouse gene annotation
  "airway",        # example dataset (human)
  "pasilla"        # example dataset (Drosophila)
)
optional <- c("readxl", "ashr")

installed <- function(p) requireNamespace(p, quietly = TRUE)

if (!check_only) {
  missing_cran <- cran_pkgs[!vapply(cran_pkgs, installed, logical(1))]
  if (length(missing_cran)) {
    message("Installing from CRAN: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran)
  }
  if (!installed("BiocManager")) install.packages("BiocManager")
  missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, installed, logical(1))]
  if (length(missing_bioc)) {
    message("Installing from Bioconductor ", as.character(BiocManager::version()), ": ",
            paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
  }
  if (installed("msigdbr")) {
    message("Fetching the MSigDB gene-set database (first use only) ...")
    ok <- tryCatch({ invisible(msigdbr::msigdbr(species = "Homo sapiens", collection = "H")); TRUE },
                   error = function(e) { message("  Could not fetch MSigDB: ", conditionMessage(e)); FALSE })
    if (ok) message("  MSigDB gene sets are ready.")
  }
}

all_pkgs <- c(cran_pkgs, bioc_pkgs)
status <- data.frame(
  package = all_pkgs,
  source = c(rep("CRAN", length(cran_pkgs)), rep("Bioconductor", length(bioc_pkgs))),
  version = vapply(all_pkgs, function(p) if (installed(p)) as.character(packageVersion(p)) else "MISSING", character(1)),
  required = ifelse(all_pkgs %in% optional, "optional", "required"),
  row.names = NULL
)
cat("\n", R.version.string, "\n\n", sep = "")
print(status, row.names = FALSE)
pandoc <- installed("rmarkdown") && rmarkdown::pandoc_available()
cat(sprintf("\npandoc (for the HTML report): %s\n", if (pandoc) "available" else "NOT FOUND (install pandoc, or use RStudio)"))
bad <- status$package[status$version == "MISSING" & status$required == "required"]
if (length(bad)) {
  cat("\nMissing required packages: ", paste(bad, collapse = ", "), "\n", sep = "")
  quit(save = "no", status = 1)
}
cat("\nAll required packages are installed.\n")
