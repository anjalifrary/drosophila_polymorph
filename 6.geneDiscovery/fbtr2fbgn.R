### FROM ALAN

## Convert FlyBase transcript IDs (FBtr) to gene IDs (FBgn) offline.
##
## Download once from FlyBase (Genes section of the precomputed files):
##   https://flybase.org/releases/current/precomputed_files/genes/
##   fbgn_fbtr_fbpp_expanded_fb_<release>.tsv.gz   (also gives gene symbols)
##   or fbgn_fbtr_fbpp_fb_<release>.tsv.gz         (IDs only)
##
## Requires: data.table

library(data.table)

## Read the FlyBase mapping file into a table: fbtr, fbgn (+ gene_symbol if present)
read_fbtr_map <- function(path) {
  l <- readLines(gzfile(path))
  cmt <- startsWith(l, "#")
  hdr <- tail(grep("\t", l[cmt], value = TRUE), 1)          # header = last commented line with tabs
  x <- fread(text = l[!cmt & nzchar(l)], header = FALSE, sep = "\t", colClasses = "character")
  if (length(hdr)) {
    nm <- trimws(sub("^#+", "", strsplit(hdr, "\t")[[1]]))
    if (length(nm) == ncol(x)) setnames(x, nm)
  }
  ## find ID columns by content, so header naming differences don't matter
  pick <- function(re) names(x)[which(vapply(x, function(v) mean(grepl(re, v)) > 0.9, TRUE))[1]]
  fbgn_col <- pick("^FBgn[0-9]+$"); fbtr_col <- pick("^FBtr[0-9]+$")
  if (is.na(fbgn_col) || is.na(fbtr_col)) stop("couldn't find FBgn/FBtr columns in ", path)
  out <- x[, .(fbtr = get(fbtr_col), fbgn = get(fbgn_col))]
  sym <- grep("^gene_symbol$", names(x), ignore.case = TRUE, value = TRUE)
  if (length(sym)) out[, gene_symbol := x[[sym[1]]]]
  unique(out[fbtr != ""])
}

## ids: character vector of FBtr IDs. Returns one row per input ID, in the same order.
fbtr_to_fbgn <- function(ids, map = "fbgn_fbtr_fbpp_expanded_fb_2026_03.tsv.gz") {
  if (is.character(map) && length(map) == 1) map <- read_fbtr_map(map)
  res <- map[data.table(fbtr = ids), on = "fbtr", mult = "first"]
  miss <- unique(res[is.na(fbgn), fbtr])
  if (length(miss))
    message(length(miss), " FBtr IDs not in this release (retired, or from another release): ",
            paste(head(miss, 10), collapse = ", "), if (length(miss) > 10) ", ...")
  res[]
}

## Example:
## map <- read_fbtr_map("fbgn_fbtr_fbpp_expanded_fb_2026_03.tsv.gz")   # read once
## variants[, fbgn := fbtr_to_fbgn(transcript, map)$fbgn]
## variants[, gene_symbol := fbtr_to_fbgn(transcript, map)$gene_symbol]