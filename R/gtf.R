# =============================================================================
# R/gtf.R  --  reading a GTF annotation file with base R and data.table
# =============================================================================
# Shared by later/ (U03 tximport, U05 sashimi plots) and the isoforms/ module
# (transcript-to-gene mapping). Source after R/utils.R. No Bioconductor needed.
# =============================================================================

# ---- GTF -------------------------------------------------------------------------------------
# Transcript and exon records of a GTF (plain or .gz) as a data.frame, read with base R.
read_gtf <- function(path, features = c("transcript", "exon")) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  type <- sub("^[^\t]*\t[^\t]*\t([^\t]*)\t.*$", "\\1", lines)
  lines <- lines[type %in% features]
  f <- data.table::tstrsplit(lines, "\t", fixed = TRUE)
  if (length(f) < 9) fail("The GTF file has no ", paste(features, collapse = "/"), " lines with 9 columns: ", path)
  attrs <- f[[9]]
  attr_val <- function(key) {  # vectorised: one regex over all lines
    m <- regexpr(paste0("(^|;\\s*)", key, ' "[^"]*"'), attrs)
    out <- rep(NA_character_, length(attrs))
    out[m > 0] <- sub(paste0('^.*', key, ' "([^"]*)"$'), "\\1", regmatches(attrs, m))
    out
  }
  data.frame(chr = f[[1]], feature = f[[3]], start = as.integer(f[[4]]), end = as.integer(f[[5]]), strand = f[[7]],
             gene_id = attr_val("gene_id"), gene_name = attr_val("gene_name"), transcript_id = attr_val("transcript_id"),
             stringsAsFactors = FALSE)
}

# transcript_id -> gene_id, gene_name (one row per transcript)
tx2gene_from_gtf <- function(path) {
  g <- read_gtf(path, features = c("transcript", "exon"))
  g <- g[!is.na(g$transcript_id), c("transcript_id", "gene_id", "gene_name")]
  g[!duplicated(g$transcript_id), ]
}
