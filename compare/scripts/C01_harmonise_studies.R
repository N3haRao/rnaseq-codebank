#!/usr/bin/env Rscript
# =============================================================================
# C01_harmonise_studies.R  --  read each study's results and match their genes
# =============================================================================
# For every study (a code-bank run or a DE results table) this step:
#   - reads the gene table and what was compared (test vs reference group)
#   - flips the sign of the fold changes when asked (flip: true)
#   - calls significant genes with ONE rule for all studies (padj and |log2FC|)
#   - picks how genes are matched across studies (gene ID, Entrez ID, symbol, or
#     human-mouse one-to-one orthologs) and gives every gene a shared key
#   - picks one ranking metric for all studies (Wald statistic or signed p value)
# Outputs: studies/<study>.csv (standard columns), matched_genes.csv (one row per
# shared gene), studies.json, validation_report.json, resolved_params.yaml.
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "compare/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "compare", "R", "compare_utils.R"))
suppressPackageStartupMessages(library(optparse))

STEP <- "C01_harmonise_studies"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character")),
                               description = "Comparison step C01: read and harmonise the studies."))
if (is_blank(opt[["run_dir"]])) fail("--run_dir is required.")
params <- load_compare_params(opt[["run_dir"]], use_resolved = FALSE)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
dir.create(file.path(out_dir, "studies"), showWarnings = FALSE)

checks <- list()
add_check <- function(id, status, message) {
  checks[[length(checks) + 1]] <<- list(id = id, status = status, message = message)
  f <- switch(status, fail = log_warn, warning = log_warn, log_info)
  f(sprintf("[%s] %s: %s", toupper(status), id, message))
}
finish <- function(ok, extra = list()) {
  write_json_file(c(list(status = if (ok) "pass" else "fail", checks = checks), extra),
                  file.path(out_dir, "validation_report.json"))
  write_step_provenance(run_dir, STEP, started, params_used = params[c("match_genes_by", "padj_cutoff", "log2fc_cutoff",
                                                                     "use_study_cutoffs", "rank_by")],
                        packages = c("babelgene"), outputs = list.files(out_dir))
  if (!ok) fail("The studies cannot be compared: see ", file.path(out_dir, "validation_report.json"))
}

# ---- 1. The list of studies ---------------------------------------------------------
studies <- params$studies %||% list()
if (length(studies) < 2) { add_check("n_studies", "fail", "At least two studies are needed for a comparison."); finish(FALSE) }
if (length(studies) > 8) { add_check("n_studies", "fail", sprintf("At most 8 studies can be compared at once (%d given).", length(studies))); finish(FALSE) }
labels <- vapply(seq_along(studies), function(i) as.character(studies[[i]]$label %||% paste0("study", i)), character(1))
if (anyDuplicated(safe_label(labels))) {
  add_check("study_labels", "fail", paste0("Study labels must be different from each other: ", paste(labels, collapse = ", "), "."))
  finish(FALSE)
}
add_check("n_studies", "pass", sprintf("%d studies: %s.", length(studies), paste(labels, collapse = ", ")))

# ---- 2. Read each study -----------------------------------------------------------
num <- function(x) suppressWarnings(as.numeric(x))
read_run_study <- function(s, dir) {
  de_path <- file.path(dir, "04_differential_expression", "de_results_all.csv")
  if (!file.exists(de_path)) stop("no 04_differential_expression/de_results_all.csv (is it a finished pipeline run?)")
  de <- data.table::fread(de_path, data.table = FALSE, colClasses = list(character = c("gene_id", "gene_symbol", "entrez_id")))
  rp <- yaml::read_yaml(file.path(dir, "01_validate_inputs", "resolved_params.yaml"))
  ds <- read_json_file(file.path(dir, "04_differential_expression", "de_summary.json")) %||% list()
  gsea <- list.files(file.path(dir, "06_enrichment"), pattern = "^gsea\\.csv$", recursive = TRUE)
  list(table = data.frame(gene_id = de$gene_id, gene_symbol = de$gene_symbol, entrez_id = de$entrez_id,
                          log2FoldChange = num(de$log2FoldChange), stat = num(de$stat), pvalue = num(de$pvalue),
                          padj = num(de$padj), significant_own = as.logical(de$significant), stringsAsFactors = FALSE),
       source = "run", organism = rp$organism %||% "other", gene_id_type = rp$gene_id_type %||% "other",
       comparison = ds$comparison %||% sprintf("%s vs %s", rp$test_level, rp$reference_level),
       test_level = rp$test_level, reference_level = rp$reference_level, condition = rp$condition,
       own_rule = ds$significance_rule, shrinkage = ds$lfc_shrinkage_used %||% "unknown",
       gsea_collections = sub("/gsea\\.csv$", "", gsea))
}
read_table_study <- function(s, path) {
  if (is.null(path) || !file.exists(path)) stop("de_file not found: ", path %||% "(none given)")
  df <- read_table_file(path)
  cols <- names(df)
  g_id <- find_column(cols, c("gene_id", "geneid", "gene", "id", "ensembl", "ensembl_id", "ensembl_gene_id", "entrez", "entrez_id", "V1"))
  g_sym <- find_column(cols, c("gene_symbol", "symbol", "gene_name", "genesymbol", "external_gene_name", "hgnc_symbol", "mgi_symbol"))
  if (is.null(g_id) && is.null(g_sym) && is.character(df[[1]])) g_id <- cols[1]
  lfc <- find_column(cols, c("log2FoldChange", "logFC", "log2FC", "log2_fold_change", "log2fc", "lfc"))
  padj <- find_column(cols, c("padj", "adj.P.Val", "FDR", "qvalue", "q_value", "p_adj", "padjust", "adj_pvalue", "BH"))
  pv <- find_column(cols, c("pvalue", "P.Value", "PValue", "p_value", "pval", "p"))
  st <- find_column(cols, c("stat", "wald_stat", "t", "statistic"))
  missing <- c(if (is.null(g_id) && is.null(g_sym)) "a gene column (gene_id or gene_symbol)",
               if (is.null(lfc)) "log2FoldChange (or logFC)", if (is.null(padj)) "padj (or FDR / adj.P.Val)")
  if (length(missing)) stop("the table has no ", paste(missing, collapse = ", "), ". Columns found: ", paste(head(cols, 15), collapse = ", "))
  ids <- as.character(df[[g_id %||% g_sym]])
  det <- detect_gene_ids(ids)
  organism <- s$organism %||% "auto"
  if (organism == "auto") organism <- if (det$organism_guess %in% c("human", "mouse")) det$organism_guess else "other"
  sym <- if (!is.null(g_sym)) as.character(df[[g_sym]]) else if (det$id_type == "symbol") ids else rep(NA_character_, length(ids))
  entrez <- if (det$id_type == "entrez") ids else rep(NA_character_, length(ids))
  if (all(is.na(sym)) && organism %in% c("human", "mouse") && det$id_type %in% c("ensembl", "entrez")) {
    ann <- tryCatch(annotate_genes(ids, organism, det$id_type), error = function(e) NULL)
    if (!is.null(ann)) { sym <- ann$gene_symbol; entrez <- ann$entrez_id }
  }
  list(table = data.frame(gene_id = ids, gene_symbol = sym, entrez_id = entrez, log2FoldChange = num(df[[lfc]]),
                          stat = if (!is.null(st)) num(df[[st]]) else NA_real_, pvalue = if (!is.null(pv)) num(df[[pv]]) else NA_real_,
                          padj = num(df[[padj]]), significant_own = NA, stringsAsFactors = FALSE),
       source = "table", organism = organism, gene_id_type = det$id_type,
       comparison = s$comparison %||% "(not stated)", test_level = NULL, reference_level = NULL, condition = NULL,
       own_rule = NULL, shrinkage = s$shrinkage %||% "unknown (external table)", gsea_collections = character(0),
       columns_used = list(gene = g_id %||% g_sym, symbol = g_sym, log2FoldChange = lfc, padj = padj, pvalue = pv, stat = st))
}

raw <- list()
ok <- TRUE
for (i in seq_along(studies)) {
  s <- studies[[i]]
  label <- labels[i]
  dir <- s$input_dir %||% s$run_dir
  r <- tryCatch(if (!is.null(s$de_file)) read_table_study(s, resolve_path(s$de_file))
                else if (!is.null(dir) && dir.exists(resolve_path(dir))) read_run_study(s, resolve_path(dir))
                else stop("neither a run folder nor a de_file was found (", dir %||% "nothing given", ")"),
                error = function(e) e)
  if (inherits(r, "error")) {
    add_check(paste0("read_", safe_label(label)), "fail", sprintf("Study '%s': %s", label, conditionMessage(r)))
    ok <- FALSE; next
  }
  r$label <- label
  r$file_label <- sprintf("%02d_%s", i, safe_label(label))
  r$input_dir <- dir
  r$flipped <- isTRUE(s$flip)
  if (r$flipped) {
    r$table$log2FoldChange <- -r$table$log2FoldChange
    r$table$stat <- -r$table$stat
    parts <- strsplit(r$comparison, " vs ", fixed = TRUE)[[1]]
    r$comparison <- if (length(parts) == 2) paste(parts[2], "vs", parts[1]) else paste0("reverse of ", r$comparison)
  }
  n_ok <- sum(!is.na(r$table$log2FoldChange))
  add_check(paste0("read_", safe_label(label)), if (n_ok >= 100) "pass" else "fail",
            sprintf("Study '%s' (%s): %d genes with a fold change; %s; organism %s; gene IDs: %s%s.", label,
                    if (r$source == "run") "code-bank run" else "results table", n_ok, r$comparison, r$organism,
                    r$gene_id_type, if (r$flipped) "; fold changes flipped" else ""))
  if (n_ok < 100) ok <- FALSE
  raw[[label]] <- r
}
if (!ok) finish(FALSE)

# ---- 3. How to match genes across studies -------------------------------------------------
orgs <- unique(vapply(raw, `[[`, "", "organism"))
id_types <- unique(vapply(raw, `[[`, "", "gene_id_type"))
entrez_cover <- vapply(raw, function(r) mean(!is.na(r$table$entrez_id) & r$table$entrez_id != ""), numeric(1))
method <- params$match_genes_by %||% "auto"
if (method == "auto") {
  method <- if (length(orgs) > 1) "ortholog"
            else if (length(id_types) == 1 && id_types %in% c("ensembl", "entrez", "symbol")) "gene_id"
            else if (all(entrez_cover > 0.5)) "entrez"
            else "symbol"
}
if (method == "ortholog") {
  if (!all(orgs %in% c("human", "mouse"))) {
    add_check("gene_matching", "fail", sprintf("Studies from different organisms (%s): only human and mouse can be matched (through orthologs). Set each study's organism, or compare studies of one organism.", paste(orgs, collapse = ", ")))
    finish(FALSE)
  }
  if (!requireNamespace("babelgene", quietly = TRUE)) {
    add_check("gene_matching", "fail", "Matching human and mouse genes needs the 'babelgene' package (install.packages('babelgene')).")
    finish(FALSE)
  }
} else if (length(orgs) > 1) {
  add_check("gene_matching", "warning", sprintf("Studies come from different organisms (%s) but genes are matched by %s, not by orthologs.", paste(orgs, collapse = ", "), method))
}

symbol_of <- function(t, type) ifelse(is.na(t$gene_symbol) | t$gene_symbol == "", if (type == "symbol") t$gene_id else NA_character_, t$gene_symbol)
ortholog_map <- NULL
if (method == "ortholog" && "mouse" %in% orgs) {
  mouse_syms <- unique(unlist(lapply(raw[vapply(raw, function(r) r$organism == "mouse", logical(1))],
                                     function(r) symbol_of(r$table, r$gene_id_type))))
  mouse_syms <- mouse_syms[!is.na(mouse_syms)]
  om <- suppressMessages(babelgene::orthologs(mouse_syms, species = "mouse", human = FALSE))
  om <- om[om$support_n >= (params$min_ortholog_support %||% 3), c("human_symbol", "symbol")]
  om <- unique(om)
  om <- om[!(om$symbol %in% om$symbol[duplicated(om$symbol)]) & !(om$human_symbol %in% om$human_symbol[duplicated(om$human_symbol)]), ]
  ortholog_map <- setNames(om$human_symbol, om$symbol)
  add_check("orthologs", if (length(ortholog_map) > 1000) "pass" else "warning",
            sprintf("%d mouse genes have a one-to-one human ortholog supported by at least %d databases (babelgene %s); the others are left out.",
                    length(ortholog_map), params$min_ortholog_support %||% 3, as.character(utils::packageVersion("babelgene"))))
}
key_of <- function(r) {
  t <- r$table
  switch(method,
    gene_id = if (r$gene_id_type == "ensembl") strip_ensembl_version(t$gene_id) else t$gene_id,
    entrez = ifelse(t$entrez_id == "", NA_character_, t$entrez_id),
    symbol = symbol_of(t, r$gene_id_type),
    ortholog = { s <- symbol_of(t, r$gene_id_type); if (r$organism == "mouse") unname(ortholog_map[s]) else toupper(s) })
}
method_text <- switch(method, gene_id = sprintf("gene ID (%s)", id_types[1]), entrez = "NCBI Entrez Gene ID",
                      symbol = "gene symbol", ortholog = "human gene symbol (mouse genes through one-to-one orthologs)")

# ---- 4. Ranking metric and significance --------------------------------------------------
has_stat <- vapply(raw, function(r) mean(!is.na(r$table$stat)) > 0.9, logical(1))
has_p <- vapply(raw, function(r) mean(!is.na(r$table$pvalue)) > 0.5, logical(1))
rank_by <- params$rank_by %||% "auto"
if (rank_by == "auto") rank_by <- if (all(has_stat)) "stat" else "signed_p"
if (rank_by == "stat" && !all(has_stat)) {
  add_check("ranking", "warning", sprintf("rank_by is stat, but %s %s no test statistic: signed p values are used for every study instead.",
                                          paste(names(raw)[!has_stat], collapse = ", "), if (sum(!has_stat) == 1) "has" else "have"))
  rank_by <- "signed_p"
}
rank_text <- if (rank_by == "stat") "DESeq2 Wald statistic" else "sign(log2 fold change) x -log10(p value)"
if (rank_by == "signed_p" && !all(has_p)) add_check("ranking", "warning", sprintf("%s %s no p values: their ranking uses -log10(padj) instead.",
                                                       paste(names(raw)[!has_p], collapse = ", "), if (sum(!has_p) == 1) "has" else "have"))

own <- isTRUE(params$use_study_cutoffs)
rule <- sprintf("padj < %g and |log2FC| > %g", params$padj_cutoff, params$log2fc_cutoff)
shrink <- unique(vapply(raw, `[[`, "", "shrinkage"))
if (length(shrink) > 1) add_check("fold_change_shrinkage", "warning",
  sprintf("The studies' fold changes were estimated differently (%s). Shrunken and unshrunken fold changes differ most for low-count genes, which lowers their correlation; significant-gene overlaps and RRHO are less affected.", paste(shrink, collapse = ", ")))

info <- list()
tables <- list()
for (label in names(raw)) {
  r <- raw[[label]]
  t <- r$table
  t$key <- key_of(r)
  t$gene <- ifelse(is.na(t$gene_symbol) | t$gene_symbol == "", t$gene_id, t$gene_symbol)
  t <- t[!is.na(t$key) & t$key != "" & !is.na(t$log2FoldChange), ]
  p_for_rank <- ifelse(is.na(t$pvalue), t$padj, t$pvalue)
  t$rank_metric <- if (rank_by == "stat") t$stat else sign(t$log2FoldChange) * -log10(pmax(p_for_rank, 1e-300))
  # several genes with one key (e.g. two Ensembl IDs, one symbol): keep the most significant
  t <- t[order(t$key, is.na(t$pvalue), t$pvalue, -abs(t$log2FoldChange)), ]
  n_dup <- sum(duplicated(t$key))
  t <- t[!duplicated(t$key), ]
  t$significant <- if (own && r$source == "run") !is.na(t$significant_own) & t$significant_own
                   else !is.na(t$padj) & t$padj < params$padj_cutoff & abs(t$log2FoldChange) > params$log2fc_cutoff
  t$direction <- ifelse(t$significant & t$log2FoldChange > 0, "up", ifelse(t$significant & t$log2FoldChange < 0, "down", "not significant"))
  t <- t[order(is.na(t$padj), t$padj), c("key", "gene", "gene_id", "gene_symbol", "log2FoldChange", "stat", "pvalue", "padj",
                                         "rank_metric", "significant", "direction")]
  write_table(t, file.path(out_dir, "studies", paste0(r$file_label, ".csv")))
  tables[[label]] <- t
  if (n_dup) log_info(sprintf("%s: %d genes shared a key with another gene; the most significant was kept.", label, n_dup))
  info[[label]] <- list(label = label, file_label = r$file_label, source = r$source, input_dir = r$input_dir %||% "",
                        comparison = r$comparison, organism = r$organism, gene_id_type = r$gene_id_type,
                        flipped = r$flipped, shrinkage = r$shrinkage, own_rule = r$own_rule %||% "",
                        n_genes = nrow(r$table), n_matched = nrow(t), n_duplicates_dropped = n_dup,
                        n_significant = sum(t$significant), n_up = sum(t$direction == "up"), n_down = sum(t$direction == "down"),
                        gsea_collections = paste(r$gsea_collections, collapse = ";"))
}
studies_df <- do.call(rbind, lapply(info, as.data.frame, stringsAsFactors = FALSE))
rownames(studies_df) <- NULL
add_check("gene_matching", "pass", sprintf("Genes are matched by %s; %s.", method_text,
          paste(sprintf("%s %d of %d genes", studies_df$label, studies_df$n_matched, studies_df$n_genes), collapse = ", ")))
add_check("significance_rule", "pass", if (own) "Each run's own significant genes are used (use_study_cutoffs: true); tables use padj and |log2FC| from the comparison settings."
          else sprintf("One rule for every study: %s.", rule))
for (label in names(raw)) if (info[[label]]$n_significant == 0)
  add_check(paste0("significant_", safe_label(label)), "warning", sprintf("Study '%s' has no significant genes with %s: overlaps with it will be empty; correlation and RRHO still use all genes.", label, rule))

# ---- 5. Shared genes ----------------------------------------------------------------------
pairs <- study_pairs(names(tables))
for (pr in pairs) {
  n_shared <- length(intersect(tables[[pr[1]]]$key, tables[[pr[2]]]$key))
  st <- if (n_shared < 100) "fail" else if (n_shared < 2000) "warning" else "pass"
  add_check(paste0("shared_", pair_folder(pr[1], pr[2])), st,
            sprintf("%s and %s share %d genes%s.", pr[1], pr[2], n_shared,
                    if (st == "fail") ": too few to compare (check the organism and gene IDs)" else if (st == "warning") ": few, so results are less reliable" else ""))
  if (st == "fail") ok <- FALSE
}
all_keys <- Reduce(union, lapply(tables, `[[`, "key"))
wide <- data.frame(key = all_keys, stringsAsFactors = FALSE)
for (label in names(tables)) {
  t <- tables[[label]]; m <- match(all_keys, t$key); fl <- info[[label]]$file_label
  if (!"gene" %in% names(wide)) wide$gene <- t$gene[m] else wide$gene <- ifelse(is.na(wide$gene), t$gene[m], wide$gene)
  wide[[paste0(fl, "_log2FC")]] <- t$log2FoldChange[m]
  wide[[paste0(fl, "_padj")]] <- t$padj[m]
  wide[[paste0(fl, "_direction")]] <- t$direction[m]
}
wide$n_studies_tested <- rowSums(!is.na(wide[, grep("_log2FC$", names(wide)), drop = FALSE]))
wide <- wide[order(-wide$n_studies_tested, wide$gene), ]
write_table(wide, file.path(out_dir, "matched_genes.csv"))

resolved <- list(match_method = method, match_method_text = method_text, rank_by = rank_by, rank_by_text = rank_text,
                 significance_rule = if (own) "each run's own rule" else rule, study_labels = names(tables))
write_yaml_file(resolved, file.path(out_dir, "resolved_params.yaml"), header = "Decisions made by step C01 (read by C02-C05).")
write_json_file(studies_df, file.path(out_dir, "studies_table.json"))
summary_out <- list(studies = studies_df, match_method = method_text, rank_by = rank_text, significance_rule = resolved$significance_rule,
                    n_genes_in_all_studies = sum(wide$n_studies_tested == length(tables)))
write_json_file(summary_out, file.path(out_dir, "studies.json"))
finish(ok, list(studies = studies_df, match_method = method_text, rank_by = rank_text))
log_info("Studies harmonised: ", method_text, "; ranking by ", rank_text, ".")
