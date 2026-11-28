library(data.table)
# library(stringr)
# library(AnnotationDbi)
# library(org.Dm.eg.db)

dataset <- "poolseq"

if (dataset == "poolseq") {
    shared_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/poolseq_dest/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.shared.bothMelSim.classed.rds")
    vcf <- "/scratch/ejy4bu/drosophila/DEST_remake/vcfs/filtering/mel/testGeneID/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.snpGap10.snpsOnly.repeatmasked.wmdust.ann.eff.vcf.gz"
    new_shared_dt <- "/project/berglandlab/anjali/drosophila_polymorphism/classification/poolseq_dest/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.shared.bothMelSim.classed.fbgn.rds"

} else if (dataset == "inbred") {
    shared_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.11_2026.shared.bothMelSim.0.05MAFfilter.classed.rds")
    vcf <- "/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/testGeneID/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.geneId.vcf.gz"
    new_shared_dt <- "/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.11_2026.shared.bothMelSim.0.05MAFfilter.classed.fbgn.rds"
} else {
    stop(paste0("Invalid dataset: ", dataset))
}

# safe_dt <- copy(shared_dt)

### 1. get coordinates of shared_dt
print("Get shared coordinates")
shared_coords <- unique(shared_dt[, .(
    CHROM = sub("^chr", "", chr),
    POS = as.integer(pos)
)])
setorder(shared_coords, CHROM, POS)
print(paste("Unique shared positions:", nrow(shared_coords)))


print("Getting coordinate query")
# Format regions as CHROM, START, END (1-based inclusive)
regions_text <- paste(
    sprintf(
        "%s\t%d\t%d",
        shared_coords$CHROM,
        shared_coords$POS,
        shared_coords$POS
    ),
    collapse = "\n"
)

print("Querying indexed VCF")

# Bash process substitution supplies the region list without a regions file
inner_cmd <- paste0(
    "bcftools query ",
    "-R <(printf '%s\\n' ",
    shQuote(regions_text),
    ") ",
    "-f '%CHROM\\t%POS\\t%INFO/EFF\\n' ",
    shQuote(vcf)
)

cmd <- paste("bash -c", shQuote(inner_cmd))

# Stream bcftools output directly into fread
eff_raw <- fread(cmd = cmd, sep = "\t", header = FALSE)
if (nrow(eff_raw) == 0) {
    stop("No VCF records retrieved. Check contig names, coordinates, VCF index, and bcftools.")
}
setnames(eff_raw, c("CHROM", "POS", "EFF"))


eff_raw[, POS := as.integer(POS)]

print(paste("Retrieved VCF records:", nrow(eff_raw)))


# ### 2. read vcf 
# print("Loading vcf")
# vcf_dt <- fread(vcf, skip = "#CHROM", select = c(1, 2, 8))
# print("Loaded vcf")
# setnames(vcf_dt, c("CHROM", "POS", "INFO"))

# vcf_dt[, chr_key := sub("^chr", "", CHROM)]
# vcf_dt[, POS := as.integer(POS)]

# vcf_dt <- vcf_dt[shared_coords, on = .(chr_key, POS = pos), nomatch = 0]
# vcf_dt[, chr_key := NULL]

# ### 3. Extract EFF from vcf
# vcf_dt <- vcf_dt[grepl("EFF=", INFO)]

# vcf_dt[, EFF := sub(".*EFF=", "", INFO)]
# vcf_dt[, EFF := sub(";.*", "", EFF)]

### 4. Split EFF entries
print("Splitting EFF annotations")
eff_raw <- eff_raw[!is.na(EFF) & EFF != "."]

eff_list <- strsplit(eff_raw$EFF, ",", fixed = TRUE)

eff_dt <- data.table(
    CHROM = rep(eff_raw$CHROM, lengths(eff_list)),
    POS = rep(eff_raw$POS, lengths(eff_list)),
    EFF_entry = unlist(eff_list, use.names = FALSE)
)

rm(eff_raw, eff_list)
gc()

### 5. Extract FBgn IDs
print("Extracting FBgn and transcript IDs")

eff_dt[, effect := sub("\\(.*", "", EFF_entry)]

eff_fields <- strsplit(sub("^[^(]*\\(", "", eff_dt$EFF_entry), "|", fixed = TRUE)

eff_dt[, gene_id_fbgn := vapply(
    eff_fields,
    function(x) if (length(x) >= 6) x[6] else NA_character_,
    character(1)
)]
eff_dt[, transcript_id_vcf := vapply(
    eff_fields,
    function(x) if (length(x) >= 9) x[9] else NA_character_,
    character(1)
)]

rm(eff_fields)
gc()


eff_dt <- eff_dt[grepl("^FBgn[0-9]+$", gene_id_fbgn)]
eff_dt <- eff_dt[effect %in% c("synonymous_variant", "missense_variant")]
print(paste("Filtered annotation entries:", nrow(eff_dt)))

### 6. Create chr-pos-fbgn id mapping
print("Creating chr-pos-FBgn mapping")

fbgn_by_pos <- unique(eff_dt[, .(CHROM, POS, gene_id_fbgn)])

fbgn_by_pos[, n_fbgn := uniqueN(gene_id_fbgn), by = .(CHROM, POS)]
fbgn_single <- unique(fbgn_by_pos[n_fbgn == 1, .(CHROM, POS, gene_id_fbgn)])


### 7. Update shared dt with SnpEff annotated gene mel
print("Adding unique position FBgn mappings")

fbgn_single[, chr_key := sub("^chr", "", CHROM)]
shared_dt[, chr_key := sub("^chr", "", chr)]

shared_dt[, geneID_fbgn := NA_character_]
shared_dt[fbgn_single, on = .(chr_key, pos = POS), geneID_fbgn := i.gene_id_fbgn]

shared_dt[, chr_key := NULL]
fbgn_single[, chr_key := NULL]

shared_dt[!is.na(variant.id_mel), .(
    total = .N,
    mapped = sum(!is.na(geneID_fbgn)),
    unmapped = sum(is.na(geneID_fbgn)),
    mapping_rate = round(100 * sum(!is.na(geneID_fbgn)) / .N, 2)
)]

### 8. resolve ambiguous positions that weren't mapped
print("Resolving ambiguous positions using transcript IDs")
eff_dt[, chr_key := sub("^chr", "", CHROM)]

missing <- shared_dt[!is.na(variant.id_mel) & is.na(geneID_fbgn)]
missing_pos <- missing[, .(chr_key=sub("^chr", "", chr), POS = as.integer(pos), transcript_id_mel)]
missing_transcripts <- eff_dt[missing_pos, on = .(chr_key, POS), .(
        chr_key, POS,
        transcript_id_mel = i.transcript_id_mel,
        transcript_id_vcf, gene_id_fbgn)
]
missing_transcripts <- missing_transcripts[transcript_id_mel == transcript_id_vcf]
missing_transcripts <- unique(
    missing_transcripts[, .(
        chr_key,
        POS,
        transcript_id_mel,
        gene_id_fbgn
    )]
)

print("Transcript-resolved mappings:")
print(head(missing_transcripts))
shared_dt[, chr_key := sub("^chr", "", chr)]

shared_dt[missing_transcripts, on = .(chr_key, pos = POS, transcript_id_mel = transcript_id_mel
    ), geneID_fbgn := i.gene_id_fbgn]


shared_dt[, chr_key := NULL]
eff_dt[, chr_key := NULL]


print("Mapping statistics after transcript resolution")
shared_dt[!is.na(variant.id_mel),
    .(
        total = .N,
        mapped = sum(!is.na(geneID_fbgn)),
        unmapped = sum(is.na(geneID_fbgn)),
        mapping_rate = round(100 * sum(!is.na(geneID_fbgn)) / .N, 2)
    )
]

print(head(shared_dt[!is.na(variant.id_mel) & !is.na(geneID_fbgn), .(chr, pos, transcript_id_mel, gene_mel, geneID_fbgn)]))

saveRDS(shared_dt, new_shared_dt)
print(paste("Saved to:", new_shared_dt))

stop("Stopping here for debugging")




shared_dt[
    is.na(geneID_fbgn) &
    !is.na(gene_id_fbgn_fromGFF),
    geneID_fbgn := gene_id_fbgn_fromGFF
]
shared_dt[!is.na(variant.id_mel), .(
    total = .N,
    mapped = sum(!is.na(geneID_fbgn)),
    unmapped = sum(is.na(geneID_fbgn)),
    mapping_rate = round(100 * sum(!is.na(geneID_fbgn)) / .N, 4)
)]











missing_transcripts <- eff_dt[
    missing_pos,
    on = .(CHROM = chr, POS = pos),
    .(
        CHROM,
        POS,
        transcript_id = vapply(
            strsplit(sub("^[^(]*\\(", "", EFF_entry), "|", fixed = TRUE),
            function(x) if (length(x) >= 5) x[5] else NA_character_,
            character(1)
        ),
        gene_id_fbgn
    )
]





gene_map <- unique(shared_dt[
    !is.na(gene_sim) & !is.na(geneID_fbgn),
    .(gene_sim, geneID_fbgn)
])
gene_map[, .N, by = gene_sim][N > 1]
gene_map[, .N, by = geneID_fbgn][N > 1]








# shared_dt <- copy(safe_dt)

# initialize fbgn shared gene ID column
shared_dt[, geneID_shared := NA_character_]

# build map: gene_sim <-> gene_mel ? fix this
gene_pairs <- unique(rbindlist(list(
    inbred_dt[, .(gene_mel, gene_sim)],
    poolseq_dt[, .(gene_mel, gene_sim)]
)))

gene_pairs <- gene_pairs[!is.na(gene_mel) & !is.na(gene_sim)]

# count number of mel genes for every sim gene
gene_pairs[, n_mel_per_sim := uniqueN(gene_mel), by = gene_sim]
gene_pairs[, n_sim_per_mel := uniqueN(gene_sim), by = gene_mel]
table(gene_pairs$n_mel_per_sim)
table(gene_pairs$n_sim_per_mel)

# Summary
gene_pairs[, .(
    n_pairs = .N,
    n_gene_mel = uniqueN(gene_mel),
    n_gene_sim = uniqueN(gene_sim),
    n_1to1 = sum(n_mel_per_sim == 1 & n_sim_per_mel == 1),
    n_sim_to_multiple_mel = uniqueN(gene_sim[n_mel_per_sim > 1]),
    n_mel_to_multiple_sim = uniqueN(gene_mel[n_sim_per_mel > 1])
)]


# keep only unambiguous 1 gene mel for each gene sim 
gene_pairs <- gene_pairs[n_mel_per_sim == 1]

# Map gene_mel -> FBgn
gene_mel_to_fbgn <- as.data.table(AnnotationDbi::select(
        org.Dm.eg.db, keys = unique(c(inbred_dt$gene_mel, poolseq_dt$gene_mel)),
        keytype = "SYMBOL", columns = "FLYBASE"))
setnames(gene_mel_to_fbgn, c("SYMBOL", "FLYBASE"), c("gene_mel", "gene_id_fbgn"))
gene_mel_to_fbgn <- unique(gene_mel_to_fbgn[, .(gene_mel, gene_id_fbgn)])

# Remove ambiguous gene_mel -> FBgn mappings (2 fbgn ids per 1 gene_mel id)
gene_mel_to_fbgn[, n_fbgn_per_mel := uniqueN(gene_id_fbgn), by = gene_mel]

missing_mel_annotation <- gene_mel_to_fbgn[is.na(gene_id_fbgn)]
gene_mel_to_fbgn <- gene_mel_to_fbgn[n_fbgn_per_mel == 1 & !is.na(gene_id_fbgn), .(gene_mel, gene_id_fbgn)]

### Map gene_sim to gene_mel to fbgn 
gene_sim_to_fbgn <- unique(merge(gene_pairs, gene_mel_to_fbgn, by = "gene_mel")[, .(gene_sim, gene_id_fbgn)])
gene_sim_to_fbgn[, n_fbgn := uniqueN(gene_id_fbgn), by = gene_sim]
if (gene_sim_to_fbgn[n_fbgn > 1, .N] > 0) {
    stop("Some gene_sim IDs map to multiple FBgn IDs")
}
gene_sim_to_fbgn <- gene_sim_to_fbgn[, .(gene_sim, gene_id_fbgn)]

# Annotate gene_mel first
shared_dt[gene_mel_to_fbgn, on = "gene_mel", geneID_shared := i.gene_id_fbgn]
# poolseq_dt[gene_mel_to_fbgn, on = "gene_mel", geneID_shared := i.gene_id_fbgn]

# Annotate gene_sim only where gene_mel mapping was unavailable
# Stop if the two annotation routes disagree
if (shared_dt[gene_sim_to_fbgn, on = "gene_sim",
    sum(!is.na(geneID_shared) & geneID_shared != i.gene_id_fbgn)
] > 0) {
    stop("Conflicting FBgn annotations")
}

# if (poolseq_dt[gene_sim_to_fbgn, on = "gene_sim",
#     sum(!is.na(geneID_shared) & geneID_shared != i.gene_id_fbgn)
# ] > 0) {
#     stop("Conflicting FBgn annotations in PoolSeq dataset")
# }

shared_dt[gene_sim_to_fbgn, on = "gene_sim",
    geneID_shared := fifelse(is.na(geneID_shared), i.gene_id_fbgn, geneID_shared)]
# poolseq_dt[gene_sim_to_fbgn, on = "gene_sim",
#     geneID_shared := fifelse(is.na(geneID_shared), i.gene_id_fbgn, geneID_shared)]

# Summarize mapping
shared_dt[, .(
    total = .N,
    mapped = sum(!is.na(geneID_shared)),
    unmapped = sum(is.na(geneID_shared))
)]

# poolseq_dt[, .(
#     total = .N,
#     mapped = sum(!is.na(geneID_shared)),
#     unmapped = sum(is.na(geneID_shared))
# )]

# Inspect remaining unmapped rows
missing <- shared_dt[is.na(geneID_shared)]
# missing_poolseq <- poolseq_dt[is.na(geneID_shared)]

View(missing[, .(chr, pos, gene_mel, gene_sim, gene_id_fbgn_fromGFF, geneID_shared)])
# View(missing_poolseq[, .(chr, pos, gene_mel, gene_sim, gene_id_fbgn_fromGFF, geneID_shared)])

missing[, .(
    total = .N,
    has_gene_mel = sum(!is.na(gene_mel)),
    has_gene_sim = sum(!is.na(gene_sim)),
    has_gff_fbgn = sum(!is.na(gene_id_fbgn_fromGFF)),
    has_gene_mel_no_gff = sum(!is.na(gene_mel) & is.na(gene_id_fbgn_fromGFF)),
    no_gene_symbols = sum(is.na(gene_mel) & is.na(gene_sim))
)]

# missing_poolseq[, .(
#     total = .N,
#     has_gene_mel = sum(!is.na(gene_mel)),
#     has_gene_sim = sum(!is.na(gene_sim)),
#     has_gff_fbgn = sum(!is.na(gene_id_fbgn_fromGFF)),
#     has_gene_mel_no_gff = sum(!is.na(gene_mel) & is.na(gene_id_fbgn_fromGFF)),
#     no_gene_symbols = sum(is.na(gene_mel) & is.na(gene_sim))
# )]

# there are many genes still missing annotations, so using gff to fill in those genes:

mel_gff <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff"

gff <- fread(
        cmd = paste("grep -v '^#'", shQuote(mel_gff)), 
        sep="\t", header=F, 
        col.names=c(
            "chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"
    ))

# map chr names to dm6 normal names 
chr_map <- gff[ feature == "region" & grepl("genome=chromosome", attributes),
    .(chr_gff = chr, chr = sub(".*Name=([^;]+);.*", "\\1", attributes))]

gff[chr_map, chr_dm6 := i.chr, on = .(chr = chr_gff)]
unique(gff[!is.na(chr_dm6), .(chr, chr_dm6)]) # X, 2L, 2R, 3L, 3R, 4, Y
gff <- gff[feature == "CDS" & !is.na(chr_dm6), .(chr_dm6, start, end, strand, phase, feature, attributes)]

# gff[, gene_name := str_extract(attributes, "(?<=;gene=)[^;]+")]
gff <- gff[, gene_id_fbgn := str_extract(attributes, "FBgn[0-9]+")][!is.na(gene_id_fbgn)]
gff <- unique(gff[!is.na(gene_id_fbgn)])

### make map to link gff gene ids to gds gene ids 
gene_id_map <- unique(
    rbindlist(list(
        inbred_dt[, .(gene_id_fbgn, gene_id_fbgn_fromGFF)],
        poolseq_dt[, .(gene_id_fbgn, gene_id_fbgn_fromGFF)]
    ))[
        !is.na(gene_id_fbgn) &
        !is.na(gene_id_fbgn_fromGFF),
        .(gene_id_fbgn_fromGFF, gene_id_fbgn)
    ]
)
setnames(gene_id_map, "gene_id_fbgn", "gene_id_fbgn_fromGDS")

gene_id_map[, .(n_GDS = uniqueN(gene_id_fbgn_fromGDS)),
    by = gene_id_fbgn_fromGFF][n_GDS > 1]       # should be 0

gene_id_map[, .(n_GFF = uniqueN(gene_id_fbgn_fromGFF)),
    by = gene_id_fbgn_fromGDS][n_GFF > 1]               # should be 0

inbred_dt[
    gff, gene_id_fbgn_fromGFF := i.gene_id_fbgn,
    on = .(
        chr = chr_dm6,
        pos >= start,
        pos <= end
    ), mult = "first"
]

poolseq_dt[
    gff, gene_id_fbgn_fromGFF := i.gene_id_fbgn,
    on = .(
        chr = chr_dm6,
        pos >= start,
        pos <= end
    ), mult = "first"
]

inbred_dt[gene_id_map,
    on = .(gene_id_fbgn_fromGFF),
    geneID_shared := fifelse(
        is.na(geneID_shared),
        i.gene_id_fbgn_fromGDS,
        geneID_shared
    )
]

poolseq_dt[gene_id_map,
    on = .(gene_id_fbgn_fromGFF),
    geneID_shared := fifelse(
        is.na(geneID_shared),
        i.gene_id_fbgn_fromGDS,
        geneID_shared
    )
]





missing_inbred_after <- inbred_dt[is.na(geneID_shared)]
missing_poolseq_after <- poolseq_dt[is.na(geneID_shared)]

View(missing_inbred_after[, .(chr, pos, gene_mel, gene_sim, gene_id_fbgn_fromGFF, geneID_shared)])
View(missing_poolseq_after[, .(chr, pos, gene_mel, gene_sim, gene_id_fbgn_fromGFF, geneID_shared)])


missing_inbred_after[, .(
    total = .N,
    has_gene_mel = sum(!is.na(gene_mel)),
    has_gene_sim = sum(!is.na(gene_sim)),
    has_gff_fbgn = sum(!is.na(gene_id_fbgn_fromGFF)),
    has_gene_mel_no_gff = sum(!is.na(gene_mel) & is.na(gene_id_fbgn_fromGFF)),
    no_gene_symbols = sum(is.na(gene_mel) & is.na(gene_sim))
)]

missing_poolseq_after[, .(
    total = .N,
    has_gene_mel = sum(!is.na(gene_mel)),
    has_gene_sim = sum(!is.na(gene_sim)),
    has_gff_fbgn = sum(!is.na(gene_id_fbgn_fromGFF)),
    has_gene_mel_no_gff = sum(!is.na(gene_mel) & is.na(gene_id_fbgn_fromGFF)),
    no_gene_symbols = sum(is.na(gene_mel) & is.na(gene_sim))
)]










gene_sim_to_fbgn_1to1 <- unique(gene_sim_to_fbgn[!is.na(gene_id_fbgn), .(gene_sim, gene_mel, gene_id_fbgn)])

# check gene sim to fbgn is 1:1 mapping
gene_sim_to_fbgn_1to1[, n_fbgn_per_sim := uniqueN(gene_id_fbgn), by = gene_sim]
stopifnot(gene_sim_to_fbgn_1to1[n_fbgn_per_sim > 1, .N] == 0)
gene_sim_to_fbgn_1to1[, n_fbgn_per_sim := NULL]

### next: map gene_mel to fbgn 
mapped_fbgn_mel <- gene_mel_to_fbgn_keep$gene_id_fbgn[
    match(inbred_dt$gene_mel, gene_mel_to_fbgn_keep$gene_mel)
]


ambiguous_mel_annotation <- gene_mel_to_fbgn[n_fbgn_per_mel > 1]

gene_mel_to_fbgn_1to1 <- gene_mel_to_fbgn[n_fbgn_per_mel == 1 & !is.na(gene_id_fbgn)][, 
    .(gene_mel, gene_id_fbgn)]


# map gene sim to gene mel to gene fbgn 
gene_sim_to_fbgn <- merge(
    gene_sim_to_mel,
    gene_mel_to_fbgn_1to1,
    by = "gene_mel",
    all.x = TRUE
)

# summarize

gene_sim_to_fbgn[, .(
    n_gene_sim = uniqueN(gene_sim),
    n_gene_mel = uniqueN(gene_mel),
    n_fbgn = uniqueN(gene_id_fbgn),
    n_missing_fbgn = sum(is.na(gene_id_fbgn))
)]

ambiguous_sim
ambiguous_mel_annotation



# gene_col <- "gene_id_fbgn_fromGFF"

# gene_map <- AnnotationDbi::select(
#     org.Dm.eg.db,
#     keys = unique(shared_table$gene_mel),
#     keytype = "SYMBOL",
#     columns = c("SYMBOL", "FLYBASE")
# )

# gene_map <- as.data.table(gene_map)
# setnames(gene_map, c("SYMBOL", "FLYBASE"), c("gene_mel", "gene_id_fbgn"))

# gene_map <- unique(gene_map[, .(gene_mel, gene_id_fbgn)])



# the gff gene names are older than the gene names annotated in the gds - discrepancies are synonyms / obsolete naming

cand_classes <- c("A", "B", "F", "G", "O", "P", "X", "Y")
tsp_classes <- c("A", "B")
conv_classes <- c("F", "G", "O", "P", "X", "Y")

NS_tsp <- c("A")
NS_conv <- c("F", "O", "X")
NS_shared <- c("A", "C", "E", "F", "H", "J", "K", "N", "O", "Q", "S", "T", "W", "X")

gff <- fread(
        cmd = paste("grep -v '^#'", shQuote(mel_gff)), 
        sep="\t", header=F, 
        col.names=c(
            "chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"
    ))

# map chr names to dm6 normal names 
chr_map <- gff[ feature == "region" & grepl("genome=chromosome", attributes),
    .(chr_gff = chr, chr = sub(".*Name=([^;]+);.*", "\\1", attributes))]

gff[chr_map, chr_dm6 := i.chr, on = .(chr = chr_gff)]
unique(gff[!is.na(chr_dm6), .(chr, chr_dm6)]) # X, 2L, 2R, 3L, 3R, 4, Y
# cds_gff <- gff[feature == "CDS" & !is.na(chr_dm6), .(chr_dm6, start, end, strand, phase, attributes)]
cds_gff <- gff[!is.na(chr_dm6), .(chr_dm6, start, end, strand, phase, feature, attributes)]

cds_gff[, transcript_id := str_extract(attributes, "(?<=Parent=rna-)[^;,]+")]
cds_gff[, gene_name := str_extract(attributes, "(?<=;gene=)[^;]+")]
cds_gff[, gene_id_fbgn := str_extract(attributes, "FBgn[0-9]+")]

### make map to link gff gene ids to gds gene ids 
gene_id_map <- unique(
    rbindlist(list(
        inbred_dt[, .(gene_id_fbgn, gene_id_fbgn_fromGFF)],
        poolseq_dt[, .(gene_id_fbgn, gene_id_fbgn_fromGFF)]
    ))[
        !is.na(gene_id_fbgn) &
        !is.na(gene_id_fbgn_fromGFF),
        .(gene_id_fbgn_fromGFF, gene_id_fbgn)
    ]
)
setnames(gene_id_map, "gene_id_fbgn", "gene_id_fbgn_fromGDS")

gene_id_map[, .(n_GDS = uniqueN(gene_id_fbgn_fromGDS)),
    by = gene_id_fbgn_fromGFF][n_GDS > 1]       # should be 0

gene_id_map[, .(n_GFF = uniqueN(gene_id_fbgn_fromGFF)),
    by = gene_id_fbgn_fromGDS][n_GFF > 1]               # should be 0


# annotate poolseq & inbred datasets with gene ids from the gff based on CDS intervals 
poolseq_dt[
    cds_gff,
    test_dm6_geneID_fbgn := i.gene_id_fbgn,
    on = .(
        chr = chr_dm6,
        pos >= start,
        pos <= end
    )
]

inbred_dt[
    cds_gff,
    test_dm6_geneID_fbgn := i.gene_id_fbgn,
    on = .(
        chr = chr_dm6,
        pos >= start,
        pos <= end
    )
]

# use new GFF annotation to map back to original GDS gene ids 
poolseq_dt[gene_id_map,
    geneID_shared := i.gene_id_fbgn_fromGDS,
    on = .(test_dm6_geneID_fbgn = gene_id_fbgn_fromGFF)
]
inbred_dt[gene_id_map,
    geneID_shared := i.gene_id_fbgn_fromGDS,
    on = .(test_dm6_geneID_fbgn = gene_id_fbgn_fromGFF)
]


View(inbred_dt[is.na(geneID_shared), .(chr, pos, gene_mel, gene_sim, 
    gene_id_fbgn, gene_id_fbgn_fromGFF, test_dm6_geneID_fbgn, geneID_shared)])
unique(inbred_dt[is.na(geneID_shared), .(gene_mel, gene_sim, 
    gene_id_fbgn, gene_id_fbgn_fromGFF, test_dm6_geneID_fbgn, geneID_shared)])

inbred_dt[, test_dm6_geneID_fbgn := NULL]
poolseq_dt[, test_dm6_geneID_fbgn := NULL]

final_dt_save <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.shared.bothMelSim.classed.rds"
    sharedOnly_dt <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.shared.classed.rds"
    candidates_dt <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.candidatesABFGOPXY.classed.rds"
    