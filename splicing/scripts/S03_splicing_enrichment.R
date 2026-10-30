#!/usr/bin/env Rscript
# =============================================================================
# S03_splicing_enrichment.R  --  pathways among genes with significant splicing
# =============================================================================
# Over-representation analysis (fgsea::fora, hypergeometric test, BH within each
# collection), the same method and gene sets as step 06 of the main pipeline.
#   genes    = genes with at least one significant event (all types together,
#              genes whose events all gain / all lose inclusion, and per event type)
#   universe = genes with at least one event that passed the coverage filter
# No GSEA: splicing changes have no single gene-level ranking statistic (a gene
# can have events in both directions), so only gene lists are tested.
#
# Outputs (in <run_dir>/S03_splicing_enrichment/<collection>/):
#   ora_all_events.csv, ora_more_inclusion.csv, ora_less_inclusion.csv, ora_<EVENT>.csv,
#   *_dotplot.png ; enrichment_summary.json
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "splicing/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "R", "genesets.R"))
source(file.path(REPO_ROOT, "splicing", "R", "splicing_utils.R"))
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

STEP <- "S03_splicing_enrichment"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "Splicing step S03: pathway enrichment of spliced genes."))
if (is_blank(opt$run_dir)) fail("--run_dir is required.")
params <- load_splicing_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
log_info("Step S03: pathway enrichment of genes with significant splicing")
set.seed(params$seed)

finish_skipped <- function(reason) {
  log_warn("Enrichment skipped: ", reason)
  write_json_file(list(status = "skipped", reason = reason), file.path(out_dir, "enrichment_summary.json"))
  write_step_provenance(run_dir, STEP, started, params_used = params[c("run_enrichment", "organism")],
                        outputs = list.files(out_dir), notes = list(status = "skipped", reason = reason))
  quit(save = "no", status = 0)
}
if (!isTRUE(params$run_enrichment)) finish_skipped("run_enrichment is false.")
if (!requireNamespace("fgsea", quietly = TRUE)) fail("Package 'fgsea' is not installed.")

source_mode <- params$gene_set_source %||% "msigdbr"
use_msigdb <- source_mode %in% c("msigdbr", "both")
use_gmt <- source_mode %in% c("gmt", "both")
if (use_msigdb && !params$organism %in% c("human", "mouse")) {
  if (!use_gmt || !length(params$gmt_files)) {
    finish_skipped(sprintf("MSigDB gene sets are used for human and mouse only (organism is '%s'). Use gene_set_source: gmt with your own gene sets.", params$organism))
  }
  use_msigdb <- FALSE
}

ev <- data.table::fread(file.path(run_dir, "S02_splicing_results", "splicing_events_all.csv"), data.table = FALSE,
                        colClasses = list(character = c("gene_id", "gene_symbol")))
ev <- ev[!is.na(ev$gene_symbol) & ev$gene_symbol != "", ]
universe <- unique(ev$gene_symbol[ev$passes_coverage])
sig_ev <- ev[ev$significant, ]
dirs <- tapply(sig_ev$direction, sig_ev$gene_symbol, function(d) paste(sort(unique(d)), collapse = "+"))
gene_lists <- c(list(all_events = unique(sig_ev$gene_symbol),
                     more_inclusion = names(dirs)[dirs == "up"],      # genes whose significant events all gain inclusion in group 1
                     less_inclusion = names(dirs)[dirs == "down"]),   # ... all lose inclusion
                setNames(lapply(params$event_types, function(t) unique(ev$gene_symbol[ev$significant & ev$event_type == t])),
                         params$event_types))
if (!length(universe)) finish_skipped("No events with gene symbols passed the coverage filter.")
log_info(sprintf("Universe: %d genes with tested events; %d genes with significant events.",
                 length(universe), length(gene_lists$all_events)))

gs <- load_gene_set_collections(use_msigdb, use_gmt, params$gene_set_collections, params$gmt_files, params$organism)

run_ora <- function(sets, genes) {
  genes <- intersect(genes, universe)
  if (length(genes) < params$min_genes_for_ora) return(NULL)
  r <- as.data.frame(fgsea::fora(pathways = sets, genes = genes, universe = universe,
                                 minSize = params$gene_set_min_size, maxSize = params$gene_set_max_size))
  if (!nrow(r)) return(r)
  r$fold_enrichment <- round((r$overlap / length(genes)) / (r$size / length(universe)), 3)
  r$overlapGenes <- collapse_list(r$overlapGenes)
  r <- r[order(r$padj, r$pval), c("pathway", "pval", "padj", "overlap", "size", "fold_enrichment", "overlapGenes")]
  names(r) <- c("pathway", "pvalue", "padj", "n_overlap", "set_size", "fold_enrichment", "overlap_genes")
  r
}

plot_ora <- function(r, title, path) {
  r <- head(r[r$padj < params$enrichment_padj_cutoff, ], 15)
  if (!nrow(r)) return(FALSE)
  r$label <- factor(pretty_path(r$pathway), levels = rev(unique(pretty_path(r$pathway))))
  p <- ggplot(r, aes(fold_enrichment, label, size = n_overlap, colour = -log10(padj))) +
    geom_point() + scale_colour_gradient(low = "#9ECAE1", high = "#08306B", name = "-log10(padj)") +
    scale_size_continuous(name = "Genes") + labs(title = title, x = "Fold enrichment", y = NULL) + theme_codebank(11)
  save_plot(p, path, width = 8.5, height = max(4, 0.3 * nrow(r) + 2))
  TRUE
}

results <- list()
for (key in names(gs$collections)) {
  col <- gs$collections[[key]]
  cdir <- step_dir(out_dir, make.names(gsub(":", "_", key)))
  res <- list(source = col$source, folder = basename(cdir))
  for (lst in names(gene_lists)) {
    r <- run_ora(col$sets, gene_lists[[lst]])
    if (is.null(r)) {
      res[[lst]] <- sprintf("skipped: fewer than %d genes", params$min_genes_for_ora)
      next
    }
    write_table(r, file.path(cdir, sprintf("ora_%s.csv", lst)))
    what <- switch(lst, all_events = "all significant events",
                   more_inclusion = paste("only more inclusion in", params$group_1_name),
                   less_inclusion = paste("only less inclusion in", params$group_1_name), EVENT_LABELS[[lst]])
    plot_ora(r, sprintf("%s: genes with %s", key, tolower(what)), file.path(cdir, sprintf("ora_%s_dotplot.png", lst)))
    top <- head(r[r$padj < params$enrichment_padj_cutoff, ], 5)
    res[[lst]] <- list(n_genes = length(intersect(gene_lists[[lst]], universe)),
                       n_significant = sum(r$padj < params$enrichment_padj_cutoff, na.rm = TRUE),
                       top = if (nrow(top)) lapply(seq_len(nrow(top)), function(i) list(pathway = top$pathway[i], padj = signif(top$padj[i], 3),
                                                                                      n_overlap = top$n_overlap[i])) else list())
  }
  results[[key]] <- res
  log_info(sprintf("%s: %s pathways significant for all spliced genes.", key,
                   if (is.list(res$all_events)) res$all_events$n_significant else "skipped"))
}

write_json_file(list(
  status = "done",
  methods = list(ora = "fgsea::fora hypergeometric test; genes = genes with >= 1 significant event; universe = genes with >= 1 event passing the coverage filter",
                 p_adjust = "Benjamini-Hochberg within each collection and gene list",
                 gene_set_size = sprintf("%d to %d genes after matching", params$gene_set_min_size, params$gene_set_max_size),
                 padj_cutoff = params$enrichment_padj_cutoff,
                 no_gsea = "Splicing changes have no gene-level ranking statistic, so only over-representation is tested."),
  package_versions = pkg_versions(c("fgsea", "msigdbr", "msigdbdf")),
  universe_size = length(universe), n_genes = lapply(gene_lists, length),
  collections = results, collections_not_loaded = gs$errors
), file.path(out_dir, "enrichment_summary.json"))
write_step_provenance(run_dir, STEP, started,
                      params_used = params[c("gene_set_source", "gene_set_collections", "gmt_files", "gene_set_min_size",
                                             "gene_set_max_size", "enrichment_padj_cutoff", "min_genes_for_ora", "organism")],
                      packages = c("fgsea", "msigdbr", "ggplot2"),
                      inputs = file.path(run_dir, "S02_splicing_results", "splicing_events_all.csv"),
                      outputs = list.files(out_dir, recursive = TRUE))
log_info("Step S03 done.")
