#!/usr/bin/env Rscript
# =============================================================================
# 06_enrichment.R  --  pathway / gene-set enrichment
# =============================================================================
# What it does
#   For each gene-set collection (MSigDB via the msigdbr package, and/or your
#   own .gmt files) it runs two complementary analyses with the fgsea package:
#
#   ORA  (over-representation, fgsea::fora, hypergeometric test)
#        Are the significant UP genes, and separately the significant DOWN
#        genes, enriched in a gene set? Background (universe) = all genes that
#        were tested in step 04 and have a gene symbol.
#
#   GSEA (fgsea::fgsea, pre-ranked)
#        Uses every tested gene ranked by the DESeq2 Wald statistic, so it needs
#        no significance cutoff. NES > 0: the set is shifted towards higher
#        expression in the test group; NES < 0: towards the reference group.
#
#   Gene sets are matched by gene symbol. Sets with fewer than gene_set_min_size
#   or more than gene_set_max_size tested genes are skipped. P-values are
#   Benjamini-Hochberg adjusted within each collection and analysis.
#
# MSigDB collections (gene_set_collections)
#   H        Hallmark gene sets
#   GO:BP    Gene Ontology biological process (MSigDB C5)
#   GO:MF, GO:CC   other GO branches (MSigDB C5)
#   REACTOME Reactome pathways (MSigDB C2 CP:REACTOME)
#   KEGG     KEGG pathways (MSigDB C2 CP:KEGG_LEGACY, or CP:KEGG in older msigdbr)
#   WIKIPATHWAYS  WikiPathways (MSigDB C2 CP:WIKIPATHWAYS)
#   Mouse uses the human MSigDB with msigdbr's ortholog mapping, so collection
#   names are the same for both organisms.
#
# Outputs (in <run_dir>/06_enrichment/<collection>/)
#   ora_up.csv, ora_down.csv, gsea.csv
#   ora_up_dotplot.png, ora_down_dotplot.png, gsea_barplot.png
# and <run_dir>/06_enrichment/enrichment_summary.json (top hits per collection)
#
# Usage
#   Rscript scripts/06_enrichment.R --run_dir runs/my_run
#   Rscript scripts/06_enrichment.R --run_dir runs/my_run --gene_set_collections "H,REACTOME"
#   Rscript scripts/06_enrichment.R --run_dir runs/my_run --gene_set_source gmt --gmt_files my_sets.gmt
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "scripts/x"), ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "R", "genesets.R"))
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

STEP <- "06_enrichment"
started <- Sys.time()

option_list <- list(
  make_option("--run_dir", type = "character", help = "Run folder (required)"),
  make_option("--config", type = "character", help = "Parameter file (default: <run_dir>/params.yaml)"),
  make_option("--set", type = "character", help = "Extra overrides"),
  make_option("--gene_set_source", type = "character", help = "msigdbr | gmt | both"),
  make_option("--gene_set_collections", type = "character", help = "Comma-separated, e.g. H,GO:BP,REACTOME,KEGG"),
  make_option("--gmt_files", type = "character", help = "Comma-separated .gmt files"),
  make_option("--enrichment_padj_cutoff", type = "double", help = "Pathway padj threshold (default 0.05)")
)
opt <- parse_args(OptionParser(option_list = option_list, description = "Step 06: pathway enrichment."))
params <- load_params(opt)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
log_info("Step 06: pathway enrichment")
set.seed(params$seed)

finish_skipped <- function(reason) {
  log_warn("Enrichment skipped: ", reason)
  write_json_file(list(status = "skipped", reason = reason), file.path(out_dir, "enrichment_summary.json"))
  write_step_provenance(run_dir, STEP, started, params_used = params[c("run_enrichment", "organism")],
                        outputs = list.files(out_dir), notes = list(status = "skipped", reason = reason))
  quit(save = "no", status = 0)
}
if (!isTRUE(params$run_enrichment)) finish_skipped("run_enrichment is false.")
if (!requireNamespace("fgsea", quietly = TRUE)) fail("Package 'fgsea' is not installed (BiocManager::install('fgsea')).")

source_mode <- params$gene_set_source %||% "msigdbr"
use_msigdb <- source_mode %in% c("msigdbr", "both")
use_gmt <- source_mode %in% c("gmt", "both")
if (use_msigdb && !params$organism %in% c("human", "mouse")) {
  if (!use_gmt || !length(params$gmt_files)) {
    finish_skipped(sprintf("MSigDB gene sets are used for human and mouse only (organism is '%s'). Provide your own gene sets with gene_set_source: gmt.", params$organism))
  }
  log_warn("Organism is '", params$organism, "': skipping MSigDB, using the .gmt files only.")
  use_msigdb <- FALSE
}

# ---- DE results ---------------------------------------------------------------
de <- data.table::fread(file.path(run_dir, "04_differential_expression", "de_results_all.csv"),
                        data.table = FALSE, colClasses = list(character = c("gene_id", "gene_symbol", "entrez_id")))
symbol_col <- de$gene_symbol
if (all(is.na(symbol_col) | symbol_col == "") && params$gene_id_type == "symbol") symbol_col <- de$gene_id
de$sym <- symbol_col
tested <- de[!is.na(de$pvalue) & !is.na(de$sym) & de$sym != "", ]
if (!nrow(tested)) finish_skipped("No tested genes have gene symbols, so gene sets cannot be matched.")
universe <- unique(tested$sym)
up_genes <- unique(tested$sym[tested$direction == "up"])
down_genes <- unique(tested$sym[tested$direction == "down"])
# GSEA ranks: Wald statistic; if a symbol appears twice keep the stronger one
rk <- tested[!is.na(tested$stat), c("sym", "stat")]
rk <- rk[order(-abs(rk$stat)), ]
rk <- rk[!duplicated(rk$sym), ]
ranks <- setNames(rk$stat, rk$sym)
log_info(sprintf("Universe: %d genes with symbols; %d up, %d down; %d genes ranked for GSEA.",
                 length(universe), length(up_genes), length(down_genes), length(ranks)))

# ---- Gene sets -----------------------------------------------------------------
gs <- load_gene_set_collections(use_msigdb, use_gmt, params$gene_set_collections, params$gmt_files, params$organism)
collections <- gs$collections
load_errors <- gs$errors

# ---- Analysis helpers -----------------------------------------------------------
run_ora <- function(sets, genes) {
  if (length(genes) < params$min_genes_for_ora) return(NULL)
  r <- fgsea::fora(pathways = sets, genes = genes, universe = universe,
                   minSize = params$gene_set_min_size, maxSize = params$gene_set_max_size)
  r <- as.data.frame(r)
  if (!nrow(r)) return(r)
  r$fold_enrichment <- (r$overlap / length(intersect(genes, universe))) / (r$size / length(universe))
  r$overlapGenes <- collapse_list(r$overlapGenes)
  r <- r[order(r$padj, r$pval), c("pathway", "pval", "padj", "overlap", "size", "fold_enrichment", "overlapGenes")]
  names(r) <- c("pathway", "pvalue", "padj", "n_overlap", "set_size", "fold_enrichment", "overlap_genes")
  r$fold_enrichment <- round(r$fold_enrichment, 3)
  r
}

run_gsea <- function(sets) {
  r <- suppressWarnings(fgsea::fgsea(pathways = sets, stats = ranks,
                                     minSize = params$gene_set_min_size, maxSize = params$gene_set_max_size,
                                     eps = 0))
  r <- as.data.frame(r)
  if (!nrow(r)) return(r)
  r$leadingEdge <- collapse_list(r$leadingEdge)
  r$direction <- ifelse(r$NES > 0, paste("higher in", params$test_level), paste("higher in", params$reference_level))
  r <- r[order(r$padj, r$pval), c("pathway", "pval", "padj", "ES", "NES", "size", "direction", "leadingEdge")]
  names(r) <- c("pathway", "pvalue", "padj", "ES", "NES", "set_size", "direction", "leading_edge_genes")
  r$ES <- round(r$ES, 4); r$NES <- round(r$NES, 4)
  r
}

plot_ora <- function(r, title, path) {
  if (is.null(r) || !nrow(r)) return(FALSE)
  r <- head(r[r$padj < params$enrichment_padj_cutoff, ], 15)
  if (!nrow(r)) return(FALSE)
  r$label <- factor(pretty_path(r$pathway), levels = rev(unique(pretty_path(r$pathway))))
  p <- ggplot(r, aes(fold_enrichment, label, size = n_overlap, colour = -log10(padj))) +
    geom_point() + scale_colour_gradient(low = "#9ECAE1", high = "#08306B", name = "-log10(padj)") +
    scale_size_continuous(name = "Genes") +
    labs(title = title, x = "Fold enrichment", y = NULL) + theme_codebank(11)
  save_plot(p, path, width = 8.5, height = max(4, 0.3 * nrow(r) + 2))
  TRUE
}

plot_gsea <- function(r, title, path) {
  if (is.null(r) || !nrow(r)) return(FALSE)
  s <- r[r$padj < params$enrichment_padj_cutoff, ]
  if (!nrow(s)) return(FALSE)
  s <- rbind(head(s[s$NES > 0, ], 10), head(s[s$NES < 0, ], 10))
  s$label <- pretty_path(s$pathway)
  s$label <- factor(s$label, levels = unique(s$label[order(s$NES)]))
  s$dir <- ifelse(s$NES > 0, paste("Higher in", params$test_level), paste("Higher in", params$reference_level))
  p <- ggplot(s, aes(NES, label, fill = dir)) + geom_col() +
    geom_vline(xintercept = 0, colour = "grey30") +
    scale_fill_manual(values = setNames(c("#D55E00", "#0072B2"),
                                        c(paste("Higher in", params$test_level), paste("Higher in", params$reference_level))),
                      name = NULL) +
    labs(title = title, subtitle = sprintf("Gene sets with padj < %g (top 10 each direction)", params$enrichment_padj_cutoff),
         x = "Normalised enrichment score (NES)", y = NULL) +
    theme_codebank(11) + theme(legend.position = "bottom")
  save_plot(p, path, width = 8.5, height = max(4, 0.3 * nrow(s) + 2.5))
  TRUE
}

top_hits <- function(r, n = 5, cols) {
  if (is.null(r) || !nrow(r)) return(list())
  r <- head(r[r$padj < params$enrichment_padj_cutoff, ], n)
  if (!nrow(r)) return(list())
  lapply(seq_len(nrow(r)), function(i) {
    x <- as.list(r[i, cols]); x$padj <- signif(x$padj, 3); x
  })
}

# ---- Run ------------------------------------------------------------------------
results_summary <- list()
for (key in names(collections)) {
  col <- collections[[key]]
  cdir <- step_dir(out_dir, make.names(gsub(":", "_", key)))
  matched <- vapply(col$sets, function(s) sum(s %in% universe), integer(1))
  n_eligible <- sum(matched >= params$gene_set_min_size & matched <= params$gene_set_max_size)
  log_info(sprintf("%s: %d gene sets, %d within size limits after matching.", key, length(col$sets), n_eligible))

  ora_up <- run_ora(col$sets, up_genes)
  ora_down <- run_ora(col$sets, down_genes)
  gsea <- run_gsea(col$sets)
  if (!is.null(ora_up)) write_table(ora_up, file.path(cdir, "ora_up.csv"))
  if (!is.null(ora_down)) write_table(ora_down, file.path(cdir, "ora_down.csv"))
  write_table(gsea, file.path(cdir, "gsea.csv"))

  if (!is.null(ora_up)) plot_ora(ora_up, sprintf("%s: genes up in %s", key, params$test_level), file.path(cdir, "ora_up_dotplot.png"))
  if (!is.null(ora_down)) plot_ora(ora_down, sprintf("%s: genes down in %s", key, params$test_level), file.path(cdir, "ora_down_dotplot.png"))
  plot_gsea(gsea, sprintf("%s: GSEA, %s vs %s", key, params$test_level, params$reference_level), file.path(cdir, "gsea_barplot.png"))

  n_sig <- function(r) if (is.null(r) || !nrow(r)) 0L else sum(r$padj < params$enrichment_padj_cutoff, na.rm = TRUE)
  results_summary[[key]] <- list(
    source = col$source,
    folder = basename(cdir),
    n_gene_sets = length(col$sets),
    n_gene_sets_tested = n_eligible,
    ora_up = if (is.null(ora_up)) sprintf("skipped: fewer than %d up genes", params$min_genes_for_ora)
             else list(n_significant = n_sig(ora_up), top = top_hits(ora_up, cols = c("pathway", "padj", "n_overlap", "fold_enrichment"))),
    ora_down = if (is.null(ora_down)) sprintf("skipped: fewer than %d down genes", params$min_genes_for_ora)
               else list(n_significant = n_sig(ora_down), top = top_hits(ora_down, cols = c("pathway", "padj", "n_overlap", "fold_enrichment"))),
    gsea = list(n_significant = n_sig(gsea),
                top_higher_in_test = top_hits(gsea[gsea$NES > 0, ], cols = c("pathway", "padj", "NES")),
                top_higher_in_reference = top_hits(gsea[gsea$NES < 0, ], cols = c("pathway", "padj", "NES")))
  )
  log_info(sprintf("%s: ORA up %s, ORA down %s, GSEA %d significant.", key,
                   if (is.null(ora_up)) "skipped" else n_sig(ora_up),
                   if (is.null(ora_down)) "skipped" else n_sig(ora_down), n_sig(gsea)))
}

gs_versions <- pkg_versions(c("fgsea", "msigdbr", "msigdbdf", "msigdbrdata"))
summary <- list(
  status = "done",
  methods = list(
    ora = "fgsea::fora hypergeometric test on significant up and down genes separately; universe = tested genes with symbols",
    gsea = "fgsea::fgsea (multilevel) on all tested genes ranked by DESeq2 Wald statistic",
    p_adjust = "Benjamini-Hochberg within each collection and analysis",
    gene_set_size = sprintf("%d to %d genes after matching", params$gene_set_min_size, params$gene_set_max_size),
    padj_cutoff = params$enrichment_padj_cutoff
  ),
  package_versions = gs_versions,
  universe_size = length(universe), n_up = length(up_genes), n_down = length(down_genes),
  collections = results_summary,
  collections_not_loaded = load_errors
)
write_json_file(summary, file.path(out_dir, "enrichment_summary.json"))

write_step_provenance(run_dir, STEP, started,
                      params_used = params[c("gene_set_source", "gene_set_collections", "gmt_files", "gene_set_min_size",
                                             "gene_set_max_size", "enrichment_padj_cutoff", "min_genes_for_ora", "seed", "organism")],
                      packages = c("fgsea", "msigdbr", "msigdbdf", "msigdbrdata", "ggplot2"),
                      inputs = c(file.path(run_dir, "04_differential_expression", "de_results_all.csv"),
                                 vapply(params$gmt_files, function(g) resolve_path(g) %||% g, character(1))),
                      outputs = list.files(out_dir, recursive = TRUE))
log_info("Step 06 done.")
