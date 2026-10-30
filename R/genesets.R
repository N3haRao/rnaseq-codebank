# =============================================================================
# R/genesets.R  --  gene-set collections shared by the enrichment steps
# =============================================================================
# Used by scripts/06_enrichment.R (differential expression) and
# splicing/scripts/S03_splicing_enrichment.R (alternative splicing), so both
# load MSigDB / GMT gene sets the same way. Source after R/utils.R.
# =============================================================================

MSIGDB_MAP <- list(
  "H"            = list(collection = "H",  subcollection = NULL),
  "GO:BP"        = list(collection = "C5", subcollection = "GO:BP"),
  "GO:MF"        = list(collection = "C5", subcollection = "GO:MF"),
  "GO:CC"        = list(collection = "C5", subcollection = "GO:CC"),
  "REACTOME"     = list(collection = "C2", subcollection = "CP:REACTOME"),
  "KEGG"         = list(collection = "C2", subcollection = c("CP:KEGG_LEGACY", "CP:KEGG", "CP:KEGG_MEDICUS")),
  "WIKIPATHWAYS" = list(collection = "C2", subcollection = "CP:WIKIPATHWAYS")
)

get_msigdb <- function(key, organism) {
  spec <- MSIGDB_MAP[[key]]
  if (is.null(spec)) stop("unknown collection '", key, "'. Known: ", paste(names(MSIGDB_MAP), collapse = ", "))
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("package 'msigdbr' is not installed")
  species <- if (organism == "mouse") "Mus musculus" else "Homo sapiens"
  new_api <- "collection" %in% names(formals(msigdbr::msigdbr))
  subs <- spec$subcollection %||% list(NULL)
  last_err <- NULL
  for (sub in subs) {
    gs <- tryCatch(suppressMessages(
      if (new_api) msigdbr::msigdbr(species = species, collection = spec$collection, subcollection = sub)
      else msigdbr::msigdbr(species = species, category = spec$collection, subcategory = sub)),
      error = function(e) { last_err <<- e; NULL })
    if (!is.null(gs) && nrow(gs) > 0) {
      sets <- split(gs$gene_symbol, gs$gs_name)
      return(list(sets = lapply(sets, unique), detail = paste(c(spec$collection, sub), collapse = " ")))
    }
  }
  stop("could not load MSigDB ", key, if (!is.null(last_err)) paste0(" (", conditionMessage(last_err), ")"))
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  sets <- lapply(parts, function(p) unique(p[-(1:2)][nzchar(p[-(1:2)])]))
  names(sets) <- vapply(parts, `[`, character(1), 1)
  sets
}

# Load the requested collections. Returns list(collections = list(<name> = list(sets, source)),
# errors = list(<name> = message)). `organism` must be human or mouse for MSigDB.
load_gene_set_collections <- function(use_msigdb, use_gmt, collection_keys, gmt_files, organism) {
  collections <- list()
  errors <- list()
  if (use_msigdb) {
    for (key in collection_keys) {
      r <- tryCatch(get_msigdb(key, organism), error = function(e) e)
      if (inherits(r, "error")) {
        log_warn("Gene sets ", key, ": ", conditionMessage(r))
        errors[[key]] <- conditionMessage(r)
      } else {
        collections[[key]] <- list(sets = r$sets, source = paste("MSigDB", r$detail))
      }
    }
  }
  if (use_gmt) {
    for (g in gmt_files) {
      gp <- resolve_path(g)
      if (!file.exists(gp)) { log_warn("GMT file not found: ", g); errors[[basename(g)]] <- "file not found"; next }
      collections[[tools::file_path_sans_ext(basename(gp))]] <- list(sets = read_gmt(gp), source = paste("GMT file", basename(gp)))
    }
  }
  if (!length(collections)) {
    msg <- paste0("No gene-set collections could be loaded",
                  if (length(errors)) paste0(" (", paste(names(errors), unlist(errors), sep = ": ", collapse = "; "), ")"), ".")
    if (use_msigdb) msg <- paste(msg, "msigdbr downloads its database on first use; check the internet connection or run scripts/00_setup_environment.R once.")
    fail(msg)
  }
  list(collections = collections, errors = errors)
}

collapse_list <- function(x) vapply(x, function(g) paste(g, collapse = ";"), character(1))

pretty_path <- function(x, width = 50) {
  x <- gsub("_", " ", sub("^(HALLMARK|GOBP|GOMF|GOCC|REACTOME|KEGG|WP|KEGG_MEDICUS)_", "", x))
  x <- tolower(x)
  x <- paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  ifelse(nchar(x) > width, paste0(substr(x, 1, width - 3), "..."), x)
}

