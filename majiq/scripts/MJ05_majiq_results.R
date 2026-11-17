#!/usr/bin/env Rscript
# =============================================================================
# MJ05_majiq_results.R  --  tables, summary and figures from voila's output
# =============================================================================
# From voila_all.tsv (one row per LSV, junction values separated by ';'):
#   per LSV the junction with the highest P(|dPSI| >= dpsi_threshold); the LSV
#   "changes" when that probability >= probability_threshold (MAJIQ's own rule).
#   Direction = sign of that junction's E[dPSI] (test - reference).
# From modulize/: number of modules per event type (cassette exon, alt 5'/3', ...).
# voila names some columns differently between versions; they are found by name
# (majiq_utils.R voila_columns), and the step stops with the list of columns it saw
# when the needed ones are missing.
#
# Outputs (in <run_dir>/MJ05_majiq_results/): majiq_lsvs_all.csv, majiq_lsvs_changing.csv,
# majiq_genes_changing.csv, modulize_event_types.csv, majiq_summary.json,
# lsv_volcano.png, event_types.png
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "majiq/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "majiq", "R", "majiq_utils.R"))
suppressPackageStartupMessages({ library(optparse); library(ggplot2) })

STEP <- "MJ05_majiq_results"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character"))))
params <- load_majiq_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
d4 <- file.path(run_dir, "MJ04_majiq_deltapsi")
ref <- params$reference_level; test <- params$test_level
P <- params$probability_threshold; Tt <- params$dpsi_threshold

t <- read_voila_tsv(file.path(d4, "voila_all.tsv"))
if (is.null(t) || !nrow(t)) fail("voila_all.tsv has no LSVs.")
cl <- voila_columns(names(t), params$ref_group, params$test_group)
if (is.na(cl$psi_ref)) cl$psi_ref <- voila_columns(names(t), ref, test)$psi_ref
if (is.na(cl$psi_test)) cl$psi_test <- voila_columns(names(t), ref, test)$psi_test
needed <- c("lsv_id", "dpsi", "prob")
if (any(is.na(unlist(cl[needed])))) {
  fail(sprintf("voila_all.tsv: could not find the column(s) for %s. Columns in the file: %s",
               paste(needed[is.na(unlist(cl[needed]))], collapse = ", "), paste(names(t), collapse = ", ")))
}
dl <- split_num(t[[cl$dpsi]]); pl <- split_num(t[[cl$prob]])
rl <- if (!is.na(cl$psi_ref)) split_num(t[[cl$psi_ref]]) else NULL
tl <- if (!is.na(cl$psi_test)) split_num(t[[cl$psi_test]]) else NULL
jl <- if (!is.na(cl$junctions)) strsplit(ifelse(is.na(t[[cl$junctions]]), "", t[[cl$junctions]]), ";") else NULL
best <- vapply(seq_len(nrow(t)), function(i) {
  p <- pl[[i]]; d <- dl[[i]]
  if (!length(p) || all(is.na(p))) return(NA_integer_)
  o <- order(-p, -abs(d)); as.integer(o[1])
}, 0L)
at <- function(l, i) if (is.null(l)) NA_real_ else vapply(seq_along(i), function(k) if (is.na(i[k]) || i[k] > length(l[[k]])) NA_real_ else l[[k]][i[k]], 0)
col_or_na <- function(k) if (is.na(cl[[k]])) rep(NA_character_, nrow(t)) else t[[cl[[k]]]]
res <- data.frame(lsv_id = t[[cl$lsv_id]], gene_id = col_or_na("gene_id"), gene_name = col_or_na("gene_name"),
                  lsv_type = col_or_na("lsv_type"), chr = col_or_na("chr"), strand = col_or_na("strand"),
                  n_junctions = lengths(dl),
                  junction = if (is.null(jl)) NA_character_ else vapply(seq_along(best), function(k) if (is.na(best[k]) || best[k] > length(jl[[k]])) NA_character_ else jl[[k]][best[k]], ""),
                  dpsi = at(dl, best), max_abs_dpsi = vapply(dl, function(d) if (all(is.na(d))) NA_real_ else max(abs(d), na.rm = TRUE), 0),
                  probability_changing = at(pl, best), psi_reference = at(rl, best), psi_test = at(tl, best),
                  stringsAsFactors = FALSE)
res$gene_name[is.na(res$gene_name) | res$gene_name == ""] <- res$gene_id[is.na(res$gene_name) | res$gene_name == ""]
res$changing <- !is.na(res$probability_changing) & res$probability_changing >= P
res$direction <- ifelse(!res$changing, "not changing", ifelse(res$dpsi > 0, "more inclusion", "less inclusion"))
res <- res[order(-res$probability_changing, -res$max_abs_dpsi, na.last = TRUE), ]
write_table(res, file.path(out_dir, "majiq_lsvs_all.csv"))
ch <- res[res$changing, ]
write_table(ch, file.path(out_dir, "majiq_lsvs_changing.csv"))
genes <- if (nrow(ch)) do.call(rbind, lapply(split(ch, ch$gene_id %||% ch$gene_name), function(d) data.frame(
  gene_id = d$gene_id[1], gene_name = d$gene_name[1], n_lsvs = nrow(d), max_abs_dpsi = max(abs(d$dpsi), na.rm = TRUE),
  max_probability = max(d$probability_changing), stringsAsFactors = FALSE))) else
  data.frame(gene_id = character(0), gene_name = character(0), n_lsvs = integer(0), max_abs_dpsi = numeric(0), max_probability = numeric(0))
genes <- genes[order(-genes$max_probability, -genes$max_abs_dpsi), , drop = FALSE]
write_table(genes, file.path(out_dir, "majiq_genes_changing.csv"))
voila_ch <- read_voila_tsv(file.path(d4, "voila_changing.tsv"))

# ---- modulize: modules per event type ------------------------------------------------------------
MOD_LABELS <- c(cassette = "cassette exon", tandem_cassette = "tandem cassette exons", alt5prime = "alternative 5' site",
                alt3prime = "alternative 3' site", alt3and5prime = "alternative 5' and 3' sites", p_alt5prime = "putative alternative 5' site",
                p_alt3prime = "putative alternative 3' site", alternative_intron = "intron retention / alternative intron",
                mutually_exclusive = "mutually exclusive exons", alternate_first_exon = "alternative first exon",
                alternate_last_exon = "alternative last exon", p_alternate_first_exon = "putative alternative first exon",
                p_alternate_last_exon = "putative alternative last exon", multi_exon_spanning = "multi-exon spanning",
                orphan_junction = "orphan junction", exitron = "exitron", other = "other")
mod_dir <- file.path(d4, "modulize")
mod <- data.frame(event_type = character(0), label = character(0), n_events = integer(0), n_changing = integer(0))
if (dir.exists(mod_dir)) {
  for (k in names(MOD_LABELS)) {
    f <- file.path(mod_dir, paste0(k, ".tsv"))
    if (!file.exists(f)) next
    x <- read_voila_tsv(f)
    n <- if (is.null(x)) 0L else nrow(x)
    chg <- grep("changing", names(x %||% list()), ignore.case = TRUE, value = TRUE)
    nc <- if (length(chg) && n) sum(toupper(x[[chg[1]]]) %in% c("TRUE", "T", "1", "YES")) else NA_integer_
    mod[nrow(mod) + 1, ] <- list(k, MOD_LABELS[[k]], n, nc)
  }
}
write_table(mod, file.path(out_dir, "modulize_event_types.csv"))

rule <- sprintf("P(|dPSI| >= %s) >= %s on at least one junction of the LSV", Tt, P)
summ <- list(
  comparison = sprintf("%s vs %s", test, ref), method = sprintf("MAJIQ %s deltapsi", params$majiq_version),
  significance_rule = rule,
  sign_convention = sprintf("dPSI = PSI(%s) - PSI(%s) (MAJIQ -psi1 %s, -psi2 %s); positive = more inclusion of that junction in %s", test, ref, ref, test, test),
  n_lsvs_quantified = nrow(res), n_lsvs_changing = nrow(ch), n_up = sum(ch$dpsi > 0, na.rm = TRUE), n_down = sum(ch$dpsi < 0, na.rm = TRUE),
  n_genes_changing = nrow(genes), n_lsvs_changing_voila = if (is.null(voila_ch)) NULL else nrow(voila_ch),
  modulize = if (nrow(mod)) "done" else "not available",
  modulize_by_type = if (nrow(mod)) mod else NULL,
  modulize_text = if (nrow(mod)) paste(sprintf("%s %d", mod$label[mod$n_events > 0], mod$n_events[mod$n_events > 0]), collapse = ", ") else NULL,
  columns_used = cl[!is.na(unlist(cl))],
  top_lsvs = head(ch[, c("gene_name", "lsv_id", "dpsi", "probability_changing")], 10)
)
write_json_file(summ, file.path(out_dir, "majiq_summary.json"))
log_info(sprintf("%d of %d LSVs change (%s) in %d genes.", nrow(ch), nrow(res), rule, nrow(genes)))

# ---- Figures ---------------------------------------------------------------------------------------
dc <- direction_colours(params)
cols <- c(`more inclusion` = dc[["up"]], `less inclusion` = dc[["down"]], `not changing` = "grey70")
v <- res[!is.na(res$dpsi) & !is.na(res$probability_changing), ]
if (nrow(v)) {
  p <- ggplot(v[order(v$changing), ], aes(dpsi, probability_changing, colour = direction)) +
    geom_point(size = 0.9, alpha = 0.7) + scale_colour_manual(values = cols) +
    geom_hline(yintercept = P, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = c(-1, 1) * Tt, linetype = "dashed", colour = "grey50") +
    labs(x = sprintf("E[dPSI] (%s - %s), junction with the highest probability", test, ref),
         y = sprintf("P(|dPSI| >= %s)", Tt), colour = NULL, title = "MAJIQ: local splicing variations") + theme_codebank()
  save_plot(p, file.path(out_dir, "lsv_volcano.png"), 7.5, 5.5)
}
if (nrow(mod) && any(mod$n_events > 0)) {
  m <- mod[mod$n_events > 0, ]
  m$label <- factor(m$label, levels = rev(m$label[order(m$n_events)]))
  p <- ggplot(m, aes(n_events, label)) + geom_col(fill = dc[["up"]], width = 0.7) +
    labs(x = "modules", y = NULL, title = "MAJIQ modulize: event types", subtitle = rule) + theme_codebank()
  save_plot(p, file.path(out_dir, "event_types.png"), 7.5, 0.6 + 0.35 * nrow(m))
}
write_step_provenance(run_dir, STEP, started, packages = "ggplot2", params_used = params[c("dpsi_threshold", "probability_threshold")],
                      outputs = list.files(out_dir))
log_info("Step MJ05 done.")
