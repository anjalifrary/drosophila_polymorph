
### Libraries
library(data.table)
library(ggplot2)
library(patchwork)



################################################################################
### 1. Prepare NLP / merge with snp table
################################################################################


load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_melanogaster.11_08_2026.nlpTable.paramask.genmap.busco.repeatmask.wmdust.Rdata")
mel_nlp <- as.data.table(nlp)
rm(nlp)


### JUST UNTIL I FIX NLP 
setnames(mel_nlp, "nLocales_poly", "nLocales_poly_proj_round")


tsp <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.shared.bothMelSim.classed.MAF.rds")
str(tsp)


################################################################################
### 2. FBtr data
################################################################################

# 1. Functions to read file (from Alan!)

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
  fbgn_col <- pick("^FBgn[0-9]+$"); fbtr_col <- pick("^FBtr[0-9]+$");
  if (is.na(fbgn_col) || is.na(fbtr_col)) stop("couldn't find FBgn/FBtr columns in ", path)
  out <- x[, .(fbtr = get(fbtr_col), fbgn = get(fbgn_col))]
  sym <- grep("^gene_symbol$", names(x), ignore.case = TRUE, value = TRUE)
  ann <- grep("^annotation_ID$", names(x), ignore.case=T, value=T)
  if (length(sym)) out[, gene_symbol := x[[sym[1]]]]
  if (length(ann)) out[, gene_annotation := x[[ann[1]]]]
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

### load map file: [ fbtr, fbgn, gene_symbol, gene_annotation ]
map <- read_fbtr_map("/scratch/ejy4bu/drosophila/fbgn_fbtr_fbpp_expanded_fb_2026_03.tsv")   # read once
mel_nlp[, fbgn := fbtr_to_fbgn(transcript, map)$fbgn]
mel_nlp[, gene_symbol := fbtr_to_fbgn(transcript, map)$gene_symbol]
mel_nlp



################################################################################
### 3. Add Enzyme Info
################################################################################

enz <- fread("/scratch/ejy4bu/drosophila/Dmel_enzyme_data_fb_2026_03.tsv", skip=4)
enz.ag <- enz[,list(enz=.N),gene_id]
mel_nlp <- merge(mel_nlp, enz.ag, by.x="fbgn", by.y="gene_id", all.x=T)


tsp.small <- tsp[, c("classification", "maf_sim", "maf_mel", "chr_dm6", "pos_dm6")]
setnames(tsp.small, c("chr_dm6", "pos_dm6"), c("chr", "pos"))
setkey(tsp.small, chr, pos)
setkey(mel_nlp, chr, pos)
mel_nlp <- merge(mel_nlp, tsp.small, all.x=T)
mel_nlp[classification%in%c("A","B"),call:="tsp"]
mel_nlp[classification%in%c("F", "G", "O", "P", "X", "Y"), call:="con"]
mel_nlp[is.na(call), call:="dmel"]


# enz counts number of enzymes associated with each gene_id in enzyme column
#   # filter for maf > 0.25
    # filter for enzymes
    # snps with nlp data, nlp >= 75
    # genmap score meaning ?? 

nlp.tsp.gene2 <- mel_nlp[maf_sim >.15 & maf_mel>.15][enz>0][!is.na(nLocales_poly_proj_round) & genmap_score==1 
# & nLocales_poly_proj_round>=75
, 
                    list(.N), 
                    list(fbgn,gene, col, call
                    # , cosmo=as.factor(nLocales_poly_proj_round>75)
                    )]

# calculates pseudo max AA position by finding max AA in nlp snps
nlp.tsp.gene.len <- mel_nlp[ , list(maxAApos = max(aa_pos)), list(gene)]

nlp.tsp.gene2 <- merge(nlp.tsp.gene2, nlp.tsp.gene.len, by="gene")

# cosmo if nlp >= 75
nlp.tsp.gene.w3 <- dcast(nlp.tsp.gene2
# [cosmo==T]
, fbgn+gene+maxAApos~col+call, value.var="N", fun.aggregate=sum)
nlp.tsp.gene.w3 <- dcast(nlp.tsp.gene2
# [cosmo==T]
, fbgn+gene~col+call, value.var="N", fun.aggregate=sum)

# this creates missense_variant[col]_tsp[call] to make column names

candidates <- nlp.tsp.gene.w3[missense_variant_tsp>=1][,c("fbgn", "gene"
# , "maxAApos"
, "missense_variant_tsp", "synonymous_variant_tsp")]$fbgn
universe <- nlp.tsp.gene.w3[missense_variant_tsp>=0][,  c("fbgn", "gene"
# , "maxAApos"
, "missense_variant_tsp", "synonymous_variant_tsp")]$fbgn
length(candidates)
length(universe)

candidates <- as.data.table(candidates)
setnames(candidates, "candidates", "gene_id")
candidates <- merge(candidates, unique(enz[, .(gene_id,gene_symbol,gene_name)]), by="gene_id", all.x=T)
candidates <- merge(candidates, nlp.tsp.gene.w3[, .(fbgn, missense_variant_tsp)], 
    by.x="gene_id", by.y="fbgn", all.x=T)


# Get coordinates of missense TSP SNPs for each gene
tsp.missense.coords <- mel_nlp[
  maf_sim > 0.15 & maf_mel > 0.15 &
    enz > 0 & !is.na(nLocales_poly_proj_round) &
    genmap_score == 1 & call == "tsp" & col == "missense_variant",
  .(tsp_missense_snps = paste(
      unique(paste(chr, pos, sep = ":")),
      collapse = ", ")), by = .(fbgn)
]

tsp.missense.flat <- mel_nlp[
  maf_sim > 0.15 &
    maf_mel > 0.15 &
    enz > 0 &
    !is.na(nLocales_poly_proj_round) &
    genmap_score == 1 &
    call == "tsp" &
    col == "missense_variant",
  unique(.SD),
  .SDcols = c("fbgn", "gene", "chr", "pos")
]

# Add coordinate lists to candidates
candidates2 <- merge(candidates, tsp.missense.coords, by.x = "gene_id", by.y = "fbgn", all.x = TRUE)

candidates3 <- merge(candidates, tsp.missense.flat, by.x = "gene_id", by.y = "fbgn", all.x = TRUE)



candidates3 <- merge(candidates3, 
    tsp[, .(chr_dm6, pos_dm6, aa_ref_mel, aa_alt_mel, aa_pos_mel, aa_pos_sim)], 
    by.x=c("chr", "pos"), by.y=c("chr_dm6", "pos_dm6"), 
    all.x=T)
setnames(candidates3, c("aa_ref_mel", "aa_alt_mel"), c("aa_ref", "aa_alt"))


################################################################################
### 4. Looking at Protein Structure
################################################################################


# install.packages("UniprotR")
library(UniprotR)

# moving to UniProt.ws

BiocManager::install("UniProt.ws")
library(UniProt.ws)

availableUniprotSpecies(pattern="Drosophila")

# for D.mel: (check taxId using above function)
up <- UniProt.ws(taxId = 7227) 


head(keytypes(up))
head(columns(up))

# figure out how to merge FBgn ids with UniProt ids (FLYBASE is a key... maybe?)
gene_ids <- unique(candidates$gene_id)

up_genes <- select(
  up,
  keys = gene_ids,
  columns = c("accession", "gene_primary", "protein_name"),
  keytype = "FLYBASE"
)


# gives: CG6415, CG4573, CG31759 - from Alan... but I only see the 3rd one in signor.
    # as nuclear encoded mito genes


# ggplot(data=mel_nlp[(!is.na(maf_mel) & !is.na(maf_sim)), ], aes(x=maf_mel, y=maf_sim)) + geom_point(size=0.2) + geom_tile() + facet_grid(~call)

# nlp.tsp.gene <- mel_nlp[!is.na(nLocales_poly_proj_round) & genmap_score==1 & busco=="Complete", 
#                     list(.N), 
#                     list(gene, col, call, cosmo=as.factor(nLocales_poly_proj_round>75))]
# nlp.tsp.gene.len <- mel_nlp[,list(maxAApos=max(aa_pos)), list(gene)]
# nlp.tsp.gene <- merge(nlp.tsp.gene, nlp.tsp.gene.len, by="gene")
