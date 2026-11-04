#!/usr/bin/env Rscript
# =============================================================================
# compare/tests/run_compare_tests.R  --  end-to-end tests of the comparison module
# =============================================================================
# Run from the repository root:  Rscript compare/tests/run_compare_tests.R
#
# Two synthetic studies are analysed with the real pipeline (DESeq2, fgsea on a
# GMT file) and a third is given as an edgeR-style results table that compared
# the groups the other way round. Planted genes, so the right answer is known:
#   genes 1-150     up in every study        genes 151-300  down in every study
#   genes 301-360   up in study 1 only       genes 361-420  up in study 1, down in study 2
# Tests:
#   A. a full run finishes and writes every expected file
#   B. genes are matched; one significance rule; the table's columns are recognised
#   C. fold changes correlate; the flipped table agrees only after flip: true (and a warning without it)
#   D. significant-gene overlaps by direction match the planted genes (Fisher tests)
#   E. RRHO: strongest overlap in the same-direction quadrants, the planted genes at the best cell;
#      rrho_map() itself on identical and on reversed rankings
#   F. pathways: the planted gene sets agree in every study (the table's GSEA is run here)
#   G. a rerun from saved parameters gives identical results; a stricter cutoff gives fewer genes
#   H. human vs mouse: matched through one-to-one orthologs
#   I. bad inputs fail with the right check (one study, missing columns, duplicate labels)
# Exit code 0 = all passed. Test runs go to compare/tests/output/ (ignored by git).
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "compare/tests/x"), "..", ".."))
})
setwd(REPO_ROOT)
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "compare", "R", "compare_utils.R"))

out_root <- file.path(REPO_ROOT, "compare", "tests", "output")
unlink(out_root, recursive = TRUE)
dir.create(out_root, recursive = TRUE)
data_dir <- file.path(out_root, "data")
dir.create(data_dir)
rscript <- file.path(R.home("bin"), "Rscript")
results <- data.frame(test = character(0), status = character(0), detail = character(0))
record <- function(test, ok, detail = "") {
  detail <- if (length(detail)) paste(detail, collapse = " ") else ""
  results[nrow(results) + 1, ] <<- list(test, if (isTRUE(ok)) "PASS" else "FAIL", detail)
  cat(sprintf("  [%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", test, if (nzchar(detail)) paste0("- ", detail) else ""))
}
run_in <- function(script, args, dir) {
  before <- list.dirs(dir, recursive = FALSE)
  log <- tempfile(fileext = ".log")
  code <- system2(rscript, c(script, "--runs_dir", shQuote(dir), args), stdout = log, stderr = log)
  after <- setdiff(list.dirs(dir, recursive = FALSE), before)
  list(code = code, run_dir = if (length(after)) after[1] else NA, log = log)
}
runs_pipe <- file.path(out_root, "pipeline_runs"); dir.create(runs_pipe)
runs_cmp <- file.path(out_root, "compare_runs"); dir.create(runs_cmp)
run_compare <- function(args) run_in("compare/scripts/run_compare.R", args, runs_cmp)
yaml_cfg <- function(x, name) { p <- file.path(data_dir, name); yaml::write_yaml(x, p); p }
js <- function(rd, f) read_json_file(file.path(rd, f))
show_log <- function(r) cat(paste0("    | ", tail(readLines(r$log, warn = FALSE), 20), "\n"), sep = "")

# ---- Synthetic studies ----------------------------------------------------------------
set.seed(11)
n_genes <- 2000
genes <- sprintf("GENE%04d", seq_len(n_genes))
up_all <- genes[1:150]; down_all <- genes[151:300]
up_1_only <- genes[301:360]; discordant <- genes[361:420]
make_study <- function(lfc, n_per_group = 4, seed) {
  set.seed(seed)
  base <- exp(rnorm(n_genes, log(300), 1))
  mu_ref <- base; mu_test <- base * 2^lfc
  disp <- 0.05
  cnt <- cbind(sapply(1:n_per_group, function(i) rnbinom(n_genes, mu = mu_ref * runif(1, 0.8, 1.2), size = 1 / disp)),
               sapply(1:n_per_group, function(i) rnbinom(n_genes, mu = mu_test * runif(1, 0.8, 1.2), size = 1 / disp)))
  colnames(cnt) <- c(paste0("ctrl", 1:n_per_group), paste0("ko", 1:n_per_group))
  list(counts = data.frame(gene_id = genes, cnt, check.names = FALSE),
       meta = data.frame(sample_id = colnames(cnt), genotype = rep(c("WT", "KO"), each = n_per_group)))
}
lfc1 <- setNames(rep(0, n_genes), genes); lfc2 <- lfc1
lfc1[up_all] <- 2; lfc1[down_all] <- -2; lfc1[up_1_only] <- 2; lfc1[discordant] <- 2
lfc2[up_all] <- 1.8; lfc2[down_all] <- -1.8; lfc2[discordant] <- -2
gmt <- file.path(data_dir, "test_sets.gmt")
other <- genes[421:n_genes]
writeLines(c(paste(c("PLANTED_UP", "test", up_all[1:60]), collapse = "\t"),
             paste(c("PLANTED_DOWN", "test", down_all[1:60]), collapse = "\t"),
             vapply(1:20, function(i) paste(c(sprintf("RANDOM_%02d", i), "test", sample(other, 40)), collapse = "\t"), "")), gmt)
pipe_cfg <- function(label, st) {
  write.csv(st$counts, file.path(data_dir, paste0(label, "_counts.csv")), row.names = FALSE)
  write.csv(st$meta, file.path(data_dir, paste0(label, "_samples.csv")), row.names = FALSE)
  yaml_cfg(list(run_name = label, counts_file = file.path(data_dir, paste0(label, "_counts.csv")),
                metadata_file = file.path(data_dir, paste0(label, "_samples.csv")), organism = "other",
                condition = "genotype", reference_level = "WT", test_level = "KO",
                gene_set_source = "gmt", gmt_files = list(gmt), gene_set_min_size = 10, lfc_shrinkage = "normal"),
           paste0(label, ".yaml"))
}
cat("\nPreparing two studies with the main pipeline\n")
r1 <- run_in("scripts/run_pipeline.R", c("--config", pipe_cfg("study1", make_study(lfc1, seed = 1))), runs_pipe)
r2 <- run_in("scripts/run_pipeline.R", c("--config", pipe_cfg("study2", make_study(lfc2, seed = 2))), runs_pipe)
record("pipeline runs for study 1 and 2", r1$code == 0 && r2$code == 0, paste(basename(c(r1$run_dir, r2$run_dir)), collapse = ", "))
if (r1$code != 0) show_log(r1)
# Study 3: an edgeR-style table that compared WT vs KO (the other way round)
de1 <- read.csv(file.path(r1$run_dir, "04_differential_expression", "de_results_all.csv"))
set.seed(3)
lfc3 <- -(lfc1 * 0.9) + rnorm(n_genes, 0, 0.25)
p3 <- ifelse(abs(lfc1) > 0, 10^-runif(n_genes, 6, 20), runif(n_genes))
tab3 <- data.frame(Gene = genes, logFC = round(lfc3, 4), logCPM = 5, PValue = p3, FDR = p.adjust(p3, "BH"))
tab3_path <- file.path(data_dir, "published_edger.tsv")
write.table(tab3, tab3_path, sep = "\t", quote = FALSE, row.names = FALSE)

studies <- list(list(label = "study1", run_dir = r1$run_dir), list(label = "study2", run_dir = r2$run_dir),
                list(label = "published", de_file = tab3_path, organism = "other", comparison = "WT vs KO", flip = TRUE))
cfg_main <- yaml_cfg(list(run_name = "test_compare", studies = studies, gmt_files = list(gmt), gene_set_source = "gmt",
                          gene_set_min_size = 10), "compare_main.yaml")

# ---- A. Full run -------------------------------------------------------------------------
cat("\nA. Full comparison of three studies\n")
a <- run_compare(c("--config", cfg_main))
rd <- a$run_dir
record("run_compare exit code 0", a$code == 0, basename(rd))
if (a$code != 0) show_log(a)
pair12 <- "study1_vs_study2"; pair13 <- "study1_vs_published"
expected <- c("params.yaml", "provenance.json", "run_summary.json",
              "C01_harmonise_studies/studies.json", "C01_harmonise_studies/validation_report.json",
              "C01_harmonise_studies/matched_genes.csv", "C01_harmonise_studies/studies/03_published.csv",
              file.path("C02_compare_genes", pair12, c("fc_scatter.png", "overlap_matrix.png", "genes_compared.csv", "overlap_table.csv")),
              "C02_compare_genes/correlation_heatmap.png", "C02_compare_genes/gene_comparison.json",
              "C02_compare_genes/genes_significant_in_several_studies.csv", "C02_compare_genes/intersections.png",
              file.path("C03_rrho", pair12, c("rrho_map.png", "rrho_grid.csv", "rrho_best_genes_up_up.csv")),
              "C03_rrho/rrho_summary.json", "C04_compare_pathways/pathway_comparison.json",
              "C04_compare_pathways/test_sets/nes_heatmap.png", file.path("C04_compare_pathways/test_sets", pair12, "nes_scatter.png"),
              "C05_compare_report/methods.md", "inputs/01_study1/04_differential_expression/de_results_all.csv")
missing <- expected[!file.exists(file.path(rd, expected))]
record("all expected outputs written", !length(missing), paste(missing, collapse = ", "))
summ <- js(rd, "run_summary.json")
record("run summary: success, all five steps", identical(summ$status, "success") && length(summ$steps) == 5)

# ---- B. Harmonising ----------------------------------------------------------------------
cat("\nB. Matching genes and one significance rule\n")
st <- js(rd, "C01_harmonise_studies/studies.json")
record("genes matched by gene ID; every study keeps its tested genes",
       grepl("gene ID", st$match_method) && all(st$studies$n_matched >= n_genes - 10), paste(st$match_method, paste(st$studies$n_matched, collapse = "/")))
record("one rule for every study", st$significance_rule == "padj < 0.05 and |log2FC| > 1")
record("edgeR column names recognised and the comparison reversed by flip",
       st$studies$comparison[3] == "KO vs WT" && isTRUE(st$studies$flipped[3]), st$studies$comparison[3])
record("ranking falls back to signed p for all (the table has no statistic)", grepl("p value", st$rank_by))

# ---- C. Fold-change correlation ------------------------------------------------------------
cat("\nC. Fold-change agreement\n")
gc <- js(rd, "C02_compare_genes/gene_comparison.json")
pr <- gc$pairs
get_pair <- function(p, a, b) p[p$study_a == a & p$study_b == b, ]
p12 <- get_pair(pr, "study1", "study2"); p13 <- get_pair(pr, "study1", "published")
record("study1 vs study2: positive Spearman correlation", p12$spearman_log2FC$r > 0.2, sprintf("%.2f", p12$spearman_log2FC$r))
record("study1 vs flipped table: strong positive correlation", p13$spearman_log2FC$r > 0.3, sprintf("%.2f", p13$spearman_log2FC$r))
noflip <- studies; noflip[[3]]$flip <- FALSE
c_nf <- run_compare(c("--config", yaml_cfg(list(run_name = "test_noflip", studies = noflip, compare_pathways = FALSE), "noflip.yaml"), "--steps", "1-2"))
gnf <- js(c_nf$run_dir, "C02_compare_genes/gene_comparison.json")
record("without flip: negative correlation and a warning to check the direction",
       get_pair(gnf$pairs, "study1", "published")$spearman_log2FC$r < -0.3 && any(grepl("flip: true", unlist(gnf$warnings))))

# ---- D. Overlaps ------------------------------------------------------------------------
cat("\nD. Overlap of significant genes\n")
gcmp <- read.csv(file.path(rd, "C02_compare_genes", pair12, "genes_compared.csv"))
uu <- gcmp$gene[gcmp$direction_study1 == "up" & gcmp$direction_study2 == "up"]
dd <- gcmp$gene[gcmp$direction_study1 == "down" & gcmp$direction_study2 == "down"]
ud <- gcmp$gene[gcmp$direction_study1 == "up" & gcmp$direction_study2 == "down"]
record("up in both: the planted genes (>= 90% found, <= 5 others)", mean(up_all %in% uu) >= 0.9 && sum(!uu %in% up_all) <= 5,
       sprintf("%d of 150, %d others", sum(up_all %in% uu), sum(!uu %in% up_all)))
record("down in both: the planted genes", mean(down_all %in% dd) >= 0.9 && sum(!dd %in% down_all) <= 5, sprintf("%d of 150", sum(down_all %in% dd)))
record("opposite directions: the planted discordant genes", mean(discordant %in% ud) >= 0.85 && p12$n_opposite_direction >= 50, sprintf("%d of 60", sum(discordant %in% ud)))
record("Fisher up/up and down/down highly significant", p12$fisher_up_up$p < 1e-50 && p12$fisher_down_down$p < 1e-50,
       sprintf("%.1g / %.1g", p12$fisher_up_up$p, p12$fisher_down_down$p))
sev <- read.csv(file.path(rd, "C02_compare_genes", "genes_significant_in_several_studies.csv"))
record("genes significant in all three studies include the planted up genes",
       mean(up_all %in% sev$gene[sev$n_up == 3]) >= 0.85, sprintf("%d", sum(up_all %in% sev$gene[sev$n_up == 3])))

# ---- E. RRHO ---------------------------------------------------------------------------
cat("\nE. Rank-rank hypergeometric overlap\n")
rr <- js(rd, "C03_rrho/rrho_summary.json")$pairs
r12 <- rr[rr$study_a == "study1" & rr$study_b == "study2", ]
record("RRHO pattern: the studies agree", grepl("agree", r12$pattern), r12$pattern)
best_uu <- read.csv(file.path(rd, "C03_rrho", pair12, "rrho_best_genes_up_up.csv"))
record("planted up genes at the strongest up/up cell", mean(up_all %in% best_uu$gene) >= 0.8, sprintf("%d of 150 in %d genes", sum(up_all %in% best_uu$gene), nrow(best_uu)))
record("discordant planted genes show up in the up/down quadrant", r12$best$up_down$neg_log10_padj > 5, sprintf("-log10 p %.1f", r12$best$up_down$neg_log10_padj))
set.seed(5); v <- setNames(rnorm(1000), sprintf("g%d", 1:1000)); w <- v + rnorm(1000, 0, 0.3)
m_same <- rrho_map(v, w); m_rev <- rrho_map(v, -w)
record("rrho_map: matching rankings -> concordant quadrants, reversed -> discordant",
       m_same$best$up_up$neg_log10_padj > 20 && m_same$best$up_down$neg_log10_padj < 1.3 && m_rev$best$up_down$neg_log10_padj > 20 && m_rev$best$up_up$neg_log10_padj < 1.3)

# ---- F. Pathways ----------------------------------------------------------------------
cat("\nF. Pathways across studies\n")
pw <- js(rd, "C04_compare_pathways/pathway_comparison.json")
record("pathways compared (runs' own GSEA + GSEA run here for the table)", identical(pw$status, "done") && "test_sets" %in% pw$collections_compared &&
         grepl("run here", pw$gsea_source$published))
allpw <- read.csv(file.path(rd, "C04_compare_pathways", "test_sets", "pathways_all_studies.csv"))
pu <- allpw[allpw$pathway == "PLANTED_UP", ]; pd <- allpw[allpw$pathway == "PLANTED_DOWN", ]
record("planted pathways significant in all three, same direction", pu$n_significant_up == 3 && pd$n_significant_down == 3)
record("random gene sets not consistently significant", sum(allpw$n_significant_up + allpw$n_significant_down >= 2 & grepl("^RANDOM", allpw$pathway)) <= 1)
meth <- readLines(file.path(rd, "C05_compare_report", "methods.md"))
record("methods text names the studies, the rule and RRHO", any(grepl("study1", meth)) && any(grepl("padj < 0.05", meth)) && any(grepl("RRHO2", meth)))

# ---- G. Reproducibility ------------------------------------------------------------------
cat("\nG. Rerun from saved parameters\n")
g1 <- run_compare(c("--based_on", rd, "--steps", "1-2"))
gg <- js(g1$run_dir, "C02_compare_genes/gene_comparison.json")
same <- function(x, y) isTRUE(all.equal(x[, c("study_a", "study_b", "n_significant_both", "up_up", "down_down")],
                                        y[, c("study_a", "study_b", "n_significant_both", "up_up", "down_down")]))
record("--based_on rerun gives identical overlaps (from its own copies of the inputs)", g1$code == 0 && same(gg$pairs, gc$pairs))
g2 <- run_compare(c("--based_on", rd, "--steps", "1-2", "--set", shQuote("padj_cutoff=0.001;log2fc_cutoff=1.5")))
g2s <- js(g2$run_dir, "C02_compare_genes/gene_comparison.json")
record("stricter cutoff: fewer genes significant in both", get_pair(g2s$pairs, "study1", "study2")$n_significant_both < p12$n_significant_both)

# ---- H. Human vs mouse -------------------------------------------------------------------
cat("\nH. Human vs mouse (orthologs)\n")
if (requireNamespace("babelgene", quietly = TRUE)) {
  od <- babelgene:::orthologs_df
  od <- od[od$taxon_id == 10090 & od$support_n >= 8 & grepl("^[A-Z][A-Z0-9]+$", od$human_symbol), c("human_symbol", "symbol")]
  od <- od[!duplicated(od$human_symbol) & !duplicated(od$symbol) &
           !(od$human_symbol %in% od$human_symbol[duplicated(od$human_symbol)]), ]
  set.seed(8); od <- od[sample(nrow(od), 1500), ]
  hl <- setNames(rnorm(1500, 0, 0.3), od$human_symbol); hl[1:100] <- 2; hl[101:200] <- -2
  ml <- unname(hl) + rnorm(1500, 0, 0.3)
  pv <- function(l) ifelse(abs(l) > 1, 1e-8, runif(length(l)))
  hp <- pv(hl); mp <- pv(ml)
  hf <- file.path(data_dir, "human.csv"); mf <- file.path(data_dir, "mouse.csv")
  write.csv(data.frame(gene_symbol = od$human_symbol, log2FoldChange = unname(hl), pvalue = hp, padj = p.adjust(hp, "BH")), hf, row.names = FALSE)
  write.csv(data.frame(symbol = od$symbol, log2FoldChange = ml, pvalue = mp, padj = p.adjust(mp, "BH")), mf, row.names = FALSE)
  hm <- run_compare(c("--config", yaml_cfg(list(run_name = "test_species", compare_pathways = FALSE,
                                                studies = list(list(label = "human", de_file = hf), list(label = "mouse", de_file = mf))),
                                           "species.yaml"), "--steps", "1-3"))
  hs <- js(hm$run_dir, "C01_harmonise_studies/studies.json")
  record("organisms detected from the symbols; matched through orthologs",
         hm$code == 0 && identical(hs$studies$organism, c("human", "mouse")) && grepl("ortholog", hs$match_method), hs$match_method)
  hp12 <- js(hm$run_dir, "C02_compare_genes/gene_comparison.json")$pairs
  record("most genes matched across species", hp12$n_shared_genes >= 1350, sprintf("%d of 1500", hp12$n_shared_genes))
  record("planted genes agree across species", hp12$up_up >= 90 && hp12$down_down >= 90 && hp12$spearman_log2FC$r > 0.5,
         sprintf("up/up %d, down/down %d, Spearman %.2f", hp12$up_up, hp12$down_down, hp12$spearman_log2FC$r))
} else record("babelgene installed (needed for human-mouse matching)", FALSE)

# ---- I. Bad inputs -------------------------------------------------------------------------
cat("\nI. Bad inputs fail clearly\n")
check_fail <- function(cfg, pattern, name) {
  r <- run_compare(c("--config", yaml_cfg(cfg, paste0(safe_label(name), ".yaml")), "--steps", "1"))
  v <- js(r$run_dir, "C01_harmonise_studies/validation_report.json")
  msg <- paste(v$checks$message[v$checks$status == "fail"], collapse = " | ")
  record(name, r$code != 0 && grepl(pattern, msg), msg)
}
check_fail(list(studies = studies[1]), "At least two studies", "one study only")
bad <- file.path(data_dir, "bad.csv"); write.csv(data.frame(gene = genes, foldchange = 1, q = 0.5), bad, row.names = FALSE)
check_fail(list(studies = list(studies[[1]], list(label = "bad", de_file = bad))), "log2FoldChange", "table without fold-change / padj columns")
check_fail(list(studies = list(studies[[1]], list(label = "study1", run_dir = r2$run_dir))), "different from each other", "duplicate labels")

cat("\n")
print(table(results$status))
write.csv(results, file.path(out_root, "test_results.csv"), row.names = FALSE)
quit(save = "no", status = if (any(results$status == "FAIL")) 1 else 0)
