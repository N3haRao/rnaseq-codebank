#!/usr/bin/env Rscript
# =============================================================================
# SU04_suppa_vs_rmats.R  --  does SUPPA2 agree with rMATS on the same samples?
# =============================================================================
# Optional (rmats_dir). Events of the two tools are matched by event type,
# chromosome, strand and their set of splice junctions ("a-b" = last base of the
# upstream exon .. first base of the downstream exon, 1-based):
#   SUPPA SE/A5/A3/MX/RI  <->  rMATS SE/A5SS/A3SS/MXE/RI   (RI: the retained intron only)
#   rMATS coordinates: *Start_0base is 0-based, *End is 1-based -> junction = (End, next Start_0base + 1)
# Sign: both dPSI = test - reference (SUPPA with the reference first; rMATS --b1 = test,
# or rmats_b1_is_test: false to flip). MX: SUPPA's inclusion form is the genomically
# first exon on both strands, rMATS's is the 1st exon on + and the 2nd on -, so SUPPA's
# dPSI is flipped for MX on the - strand. Several rows with the same junctions: the one
# with the lowest p / FDR is kept. (Ported from the lab's compare_SUPPA_rMATS.py.)
#
# Outputs (in <run_dir>/SU04_suppa_vs_rmats/): suppa_vs_rmats_matched.csv,
# significant_in_both.csv, comparison_by_type.csv, comparison_summary.json,
# suppa_vs_rmats_dpsi.png
# =============================================================================

REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  normalizePath(file.path(dirname(if (length(f)) f[1] else "suppa/scripts/x"), "..", ".."))
})
source(file.path(REPO_ROOT, "R", "utils.R"))
source(file.path(REPO_ROOT, "suppa", "R", "suppa_utils.R"))
suppressPackageStartupMessages({ library(optparse); library(ggplot2) })

STEP <- "SU04_suppa_vs_rmats"
started <- Sys.time()
opt <- parse_args(OptionParser(option_list = list(make_option("--run_dir", type = "character"))))
params <- load_suppa_params(opt$run_dir)
run_dir <- params$run_dir
out_dir <- step_dir(run_dir, STEP)
skip <- function(reason) {
  write_json_file(list(status = "skipped", reason = reason), file.path(out_dir, "comparison_summary.json"))
  log_info("Comparison with rMATS skipped: ", reason)
  quit(save = "no", status = 0)
}
rd <- resolve_path(params$rmats_dir)
if (is_blank(params$rmats_dir)) skip("no rMATS output was given (rmats_dir).")
if (is.null(rd) || !dir.exists(rd)) skip("rmats_dir was not found.")

TYPES <- c(SE = "SE", A5 = "A5SS", A3 = "A3SS", MX = "MXE", RI = "RI")
j <- function(a, b) paste0(a, "-", b)
rmats_keys <- function(rtype, r) {
  k <- switch(rtype,
    SE = cbind(j(r$upstreamEE, r$exonStart_0base + 1), j(r$exonEnd, r$downstreamES + 1)),
    A5SS = , A3SS = {
      down <- r$flankingES >= r$longExonEnd                  # flanking exon downstream of the long exon
      cbind(ifelse(down, j(r$longExonEnd, r$flankingES + 1), j(r$flankingEE, r$longExonStart_0base + 1)),
            ifelse(down, j(r$shortEE, r$flankingES + 1), j(r$flankingEE, r$shortES + 1)))
    },
    MXE = cbind(j(r$upstreamEE, r$X1stExonStart_0base + 1), j(r$X1stExonEnd, r$downstreamES + 1),
                j(r$upstreamEE, r$X2ndExonStart_0base + 1), j(r$X2ndExonEnd, r$downstreamES + 1)),
    RI = cbind(j(r$upstreamEE, r$downstreamES + 1)))
  apply(k, 1, function(x) paste(sort(x), collapse = ";"))
}
read_rmats <- function(rtype) {
  f <- file.path(rd, sprintf("%s.MATS.%s.txt", rtype, params$rmats_count_type))
  if (!file.exists(f)) return(NULL)
  r <- data.table::fread(f, data.table = FALSE)
  names(r) <- make.names(names(r))
  r <- r[!is.na(r$IncLevelDifference) & !is.na(r$FDR), , drop = FALSE]
  if (!nrow(r)) return(NULL)
  sgn <- if (isTRUE(params$rmats_b1_is_test)) 1 else -1
  t <- data.frame(chr = as.character(r$chr), strand = as.character(r$strand), key = rmats_keys(rtype, r),
                  rmats_id = r$ID, rmats_gene = gsub('"', "", r$geneSymbol), rmats_dpsi = sgn * r$IncLevelDifference,
                  rmats_fdr = r$FDR, stringsAsFactors = FALSE)
  t <- t[order(t$rmats_fdr), ]
  t[!duplicated(t[, c("chr", "strand", "key")]), ]
}

sup <- data.table::fread(file.path(run_dir, "SU03_suppa_diffsplice", "suppa_events_all.csv"), data.table = FALSE)
sup <- sup[!is.na(sup$pvalue) & !is.na(sup$dpsi) & sup$event_type %in% names(TYPES), ]
sup$key <- suppa_junction_key(sup$event_type, sup$coordinates)
sup$suppa_dpsi <- ifelse(sup$event_type == "MX" & sup$strand == "-", -sup$dpsi, sup$dpsi)
sup <- sup[order(sup$pvalue), ]
sup <- sup[!duplicated(sup[, c("event_type", "chr", "strand", "key")]), ]

cut <- params$dpsi_cutoff
rows <- list(); matched_all <- list()
for (st in names(TYPES)) {
  rt <- TYPES[[st]]
  s <- sup[sup$event_type == st, c("chr", "strand", "key", "event_id", "gene_name", "suppa_dpsi", "pvalue")]
  r <- read_rmats(rt)
  if (is.null(r)) next
  m <- merge(s, r, by = c("chr", "strand", "key"))
  m$event_type <- rt
  m$suppa_sig <- m$pvalue <= params$pvalue_cutoff & abs(m$suppa_dpsi) >= cut
  m$rmats_sig <- m$rmats_fdr <= params$rmats_fdr_cutoff & abs(m$rmats_dpsi) >= cut
  m$category <- ifelse(m$suppa_sig & m$rmats_sig, "both", ifelse(m$rmats_sig, "rMATS only", ifelse(m$suppa_sig, "SUPPA2 only", "neither")))
  m$same_direction <- sign(m$suppa_dpsi) == sign(m$rmats_dpsi)
  matched_all[[rt]] <- m
  sp <- function(d) if (nrow(d) >= 3) suppressWarnings(round(cor(d$suppa_dpsi, d$rmats_dpsi, method = "spearman"), 3)) else NA
  both <- m[m$suppa_sig & m$rmats_sig, ]
  rows[[rt]] <- data.frame(event_type = rt, suppa_tested = nrow(s), rmats_tested = nrow(r), matched = nrow(m),
                           suppa_significant = sum(s$pvalue <= params$pvalue_cutoff & abs(s$suppa_dpsi) >= cut),
                           rmats_significant = sum(r$rmats_fdr <= params$rmats_fdr_cutoff & abs(r$rmats_dpsi) >= cut),
                           significant_in_both = nrow(both), both_same_direction = sum(both$same_direction),
                           pct_rmats_significant_also_suppa = round(100 * nrow(both) / max(1, sum(m$rmats_sig)), 1),
                           spearman_dpsi_matched = sp(m),
                           spearman_dpsi_mx_minus = if (rt == "MXE") sp(m[m$strand == "-", ]) else NA)
}
if (!length(rows)) skip(sprintf("no *.MATS.%s.txt tables in rmats_dir.", params$rmats_count_type))
bt <- do.call(rbind, rows)
mall <- do.call(rbind, matched_all)
cols <- c("event_type", "gene_name", "rmats_gene", "chr", "strand", "key", "event_id", "suppa_dpsi", "pvalue", "suppa_sig",
          "rmats_id", "rmats_dpsi", "rmats_fdr", "rmats_sig", "category", "same_direction")
names(mall)[names(mall) == "key"] <- "key"
write_table(mall[, cols], file.path(out_dir, "suppa_vs_rmats_matched.csv"))
hc <- mall[mall$suppa_sig & mall$rmats_sig, cols]
write_table(hc[order(hc$rmats_fdr), ], file.path(out_dir, "significant_in_both.csv"))
write_table(bt, file.path(out_dir, "comparison_by_type.csv"))
overall_sp <- if (nrow(mall) >= 3) suppressWarnings(round(cor(mall$suppa_dpsi, mall$rmats_dpsi, method = "spearman"), 3)) else NA
write_json_file(list(status = "done", rmats_count_type = params$rmats_count_type,
                     rule = sprintf("SUPPA2 p <= %s, rMATS FDR <= %s, both |dPSI| >= %s", params$pvalue_cutoff, params$rmats_fdr_cutoff, cut),
                     n_matched = nrow(mall), n_sig_both = nrow(hc), n_sig_both_same_direction = sum(hc$same_direction),
                     n_rmats_significant_matched = sum(mall$rmats_sig), n_suppa_significant_matched = sum(mall$suppa_sig),
                     spearman_dpsi_matched = overall_sp, by_type = bt),
                file.path(out_dir, "comparison_summary.json"))
pal <- c(both = "#2a78d6", `rMATS only` = "#eb6834", `SUPPA2 only` = "#1baf7a", neither = "#c3c2b7")
p <- ggplot(mall[order(mall$category != "neither"), ], aes(rmats_dpsi, suppa_dpsi, colour = category)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 0.8, alpha = 0.7) + scale_colour_manual(values = pal) +
  facet_wrap(~event_type, nrow = 1) + coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(x = "rMATS dPSI", y = "SUPPA2 dPSI", colour = NULL,
       title = sprintf("SUPPA2 vs rMATS: %d matched events, Spearman %s", nrow(mall), overall_sp)) + theme_codebank()
save_plot(p, file.path(out_dir, "suppa_vs_rmats_dpsi.png"), 3 + 2.6 * length(unique(mall$event_type)), 4)
write_step_provenance(run_dir, STEP, started, packages = "ggplot2",
                      params_used = params[c("rmats_count_type", "rmats_fdr_cutoff", "rmats_b1_is_test", "dpsi_cutoff", "pvalue_cutoff")],
                      outputs = list.files(out_dir))
log_info(sprintf("Step SU04 done: %d events matched, %d significant in both (Spearman %s).", nrow(mall), nrow(hc), overall_sp))
