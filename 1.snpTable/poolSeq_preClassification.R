library(data.table)
library(stringr)
library(Biostrings)

#######################################
### load files
out_dir <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/shared/"

shared_dt <-readRDS(paste0(out_dir, "dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.rds"))

###############################################################################################
### 1. add AFs and NLP
###############################################################################################
temp_dt <- copy(shared_dt)

sim_nlp <-  readRDS("/scratch/ejy4bu/drosophila/DEST_remake/snpDT/nlp/sim_nlp.16_11_2026.rds")
# select cols
sim_nlp <- sim_nlp[, .(
    variant,
    global_af, poly_af,
    nLocales_poly, nLocales_total, nLocales_missing,
    nLocales_poly_proj_round, nLocales_poly_proj_binomial
)]
sim_nlp_cols <- setdiff(names(sim_nlp), "variant")
setnames(sim_nlp, sim_nlp_cols, paste0(sim_nlp_cols, "_sim"))


mel_nlp <- readRDS("/scratch/ejy4bu/drosophila/DEST_remake/snpDT/nlp/mel_nlp.18_11_2026.rds")
# select cols
mel_nlp <- mel_nlp[, .(
    variant,
    global_af, poly_af,
    nLocales_poly, nLocales_total, nLocales_missing,
    nLocales_poly_proj_round, nLocales_poly_proj_binomial
)]
mel_nlp_cols <- setdiff(names(mel_nlp), "variant")
setnames(mel_nlp, mel_nlp_cols, paste0(mel_nlp_cols, "_mel"))

# stop if duplicate variant ids 
stopifnot(uniqueN(sim_nlp$variant) == nrow(sim_nlp))
stopifnot(uniqueN(mel_nlp$variant) == nrow(mel_nlp))

temp_dt <- merge(temp_dt, sim_nlp, by.x="variant.id_sim", by.y="variant", all.x=T)

### check merge counts - SIM
uniqueN(temp_dt[!is.na(variant.id_sim), variant.id_sim]) # sim var in shared table
uniqueN(sim_nlp$variant) == uniqueN(temp_dt[!is.na(nLocales_poly_sim), variant.id_sim]) # exp: TRUE
uniqueN(temp_dt[!is.na(variant.id_sim) & is.na(nLocales_poly_sim), variant.id_sim]) # variants without nlp values
# some variants are present but no NLP because we filtered out variants if nlp == 0

temp_dt <- merge(temp_dt, mel_nlp, by.x="variant.id_mel", by.y="variant", all.x=T)

### check merge counts - MEL
uniqueN(temp_dt[!is.na(variant.id_mel), variant.id_mel]) # mel var in shared table
uniqueN(mel_nlp$variant) == uniqueN(temp_dt[!is.na(nLocales_poly_mel), variant.id_mel]) # exp: TRUE
uniqueN(temp_dt[!is.na(variant.id_mel) & is.na(nLocales_poly_mel), variant.id_mel]) # variants without nlp values
# some variants are present but no NLP because we filtered out variants if nlp == 0

stopifnot(nrow(temp_dt) == nrow(shared_dt))
shared_dt <- copy(temp_dt)
shared_dt_nlp_unfilt <- copy(temp_dt)

### filter for 5% poly_af 
temp_dt <- copy(shared_dt)
temp_dt[, poly_maf_sim := pmin(poly_af_sim, 1 - poly_af_sim)]
temp_dt[, poly_maf_mel := pmin(poly_af_mel, 1 - poly_af_mel)]

cat("\nSIM filter summary\n")
uniqueN(temp_dt[!is.na(variant.id_sim), variant.id_sim]) # number of sim variants to start
uniqueN(temp_dt[!is.na(variant.id_sim) & !is.na(poly_maf_sim), variant.id_sim]) # number of sim variants with nlp data
uniqueN(temp_dt[!is.na(variant.id_sim) & is.na(poly_maf_sim), variant.id_sim]) # number of sim variants without nlp data
uniqueN(temp_dt[!is.na(variant.id_sim) & poly_maf_sim >= 0.05, variant.id_sim]) # number of sim variants with valid MAF
uniqueN(temp_dt[!is.na(variant.id_sim) & !is.na(poly_maf_sim) & poly_maf_sim < 0.05, variant.id_sim]) # number sim variants with invalid MAF

cat("\nMEL filter summary\n")
uniqueN(temp_dt[!is.na(variant.id_mel), variant.id_mel]) 
uniqueN(temp_dt[!is.na(variant.id_mel) & !is.na(poly_maf_mel), variant.id_mel])
uniqueN(temp_dt[!is.na(variant.id_mel) & is.na(poly_maf_mel), variant.id_mel])
uniqueN(temp_dt[!is.na(variant.id_mel) & poly_maf_mel >= 0.05, variant.id_mel])
uniqueN(temp_dt[!is.na(variant.id_mel) & !is.na(poly_maf_mel) & poly_maf_mel < 0.05, variant.id_mel])


# if variant is present but missing nlp data or poly_maf < 0.05, set variant.id to NA 
temp_dt[!is.na(variant.id_sim) & (is.na(poly_maf_sim) | poly_maf_sim < 0.05), variant.id_sim := NA]
temp_dt[!is.na(variant.id_mel) & (is.na(poly_maf_mel) | poly_maf_mel < 0.05), variant.id_mel := NA]

uniqueN(temp_dt[!is.na(variant.id_sim), variant.id_sim]) # sim var retained
uniqueN(shared_dt[!is.na(variant.id_sim), variant.id_sim]) - uniqueN(temp_dt[!is.na(variant.id_sim), variant.id_sim]) # sim var removed
temp_dt <- temp_dt[!is.na(variant.id_sim) | !is.na(variant.id_mel)]

# Confirm all retained variants have valid poly-MAF
stopifnot(all(temp_dt[!is.na(variant.id_sim), poly_maf_sim] >= 0.05) & 
    all(temp_dt[!is.na(variant.id_mel), poly_maf_mel] >= 0.05))


shared_dt <- copy(temp_dt)
shared_dt_polyMAF5 <- copy(temp_dt)

### save files:
saveRDS(shared_dt_nlp_unfilt, paste0(out_dir, "dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.MAFunfiltered.rds"))
saveRDS(shared_dt_polyMAF5, paste0(out_dir, "dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.rds"))


###############################################################################################
### 2. Extract GFF info: CDS regions; chr, start, end, strand, phase 
###############################################################################################
### 2A: simulans

sim_gff <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.gff"

### gff format:
# seqid   source  feature start   end     score   strand  phase   attribute ---   
# sim_2R  Gnomon  CDS     6534278 6534771 .       -       0       ID=cds-XP_016026650.1;Parent=rna-XM_016167191.3;Dbxref=GeneID:6733679,Genbank:XP_016026650.1;Name=XP_016026650.1;gbkey=CDS;gene=LOC6733679;product=dnaJ homolog subfamily C member 13 isoform X1;protein_id=XP_016026650.1
# phase = how many bases to skip at the start of a coding sequence (CDS) feature to reach the first complete codon.

gff <- fread(
     cmd = paste("grep -v '^#'", shQuote(sim_gff)), 
     sep="\t", header=F, 
    col.names=c(
        "chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"
    ))

cds_gff <- gff[feature=="CDS", .(chr, start, end, strand, phase, attributes)]

# extract genes and transcripts
cds_gff[, transcript_id := str_extract(attributes,"(?<=Parent=rna-)[^;]+")]
cds_gff[, gene_id := str_extract(attributes, "(?<=Dbxref=GeneID:)[^,;]+")]
cds_gff[, sim_chr := chr]
cds_gff[, chr := sub("^sim_", "", sim_chr)]

sim_cds <- cds_gff[
    shared_dt[!is.na(variant.id_sim)], 
    on = .(
        sim_chr = chr_src,
        transcript_id = transcript_id_sim,
        start <= pos_src, 
        end >= pos_src
    ), nomatch = 0,
    .(# GFF information
        sim_chr = x.sim_chr,
        cds_start = x.start,
        cds_end = x.end,
        strand = x.strand,
        phase = x.phase,
        attributes = x.attributes,
        transcript_id = x.transcript_id,
        gene_id = x.gene_id,

        # shared_dt information
        chr = i.chr,
        pos = i.pos,

        variant.id_mel = i.variant.id_mel,
        ref_mel = i.ref_mel,
        alt_mel = i.alt_mel,
        effect_mel = i.effect_mel,
        codon_change_mel = i.codon_change_mel,
        aa_change_mel = i.aa_change_mel,
        aa_pos_mel = i.aa_pos_mel,

        variant.id_sim = i.variant.id_sim,
        ref_sim = i.ref_sim,
        alt_sim = i.alt_sim,
        effect_sim = i.effect_sim,
        codon_change_sim = i.codon_change_sim,
        aa_change_sim = i.aa_change_sim,
        aa_pos_sim = i.aa_pos_sim,

        # simulans coordinate/liftover information
        transcript_id_sim = i.transcript_id_sim,
        chr_src = i.chr_src,
        pos_src = i.pos_src,
        ref_src = i.ref_src,
        alt_src = i.alt_src,
        flip = i.flip,
        swap = i.swap,

        # nlp / allele frequency info
        global_af_sim = i.global_af_sim,
        poly_af_sim = i.poly_af_sim,
        nLocales_poly_sim = i.nLocales_poly_sim,
        nLocales_total_sim = i.nLocales_total_sim, 
        nLocales_missing_sim = i.nLocales_missing_sim,
        nLocales_poly_proj_round_sim = i.nLocales_poly_proj_round_sim,
        nLocales_poly_proj_binomial_sim = i.nLocales_poly_proj_binomial_sim
    )
]

sim_cds[, nt_pos_sim := as.integer(sub(".*c\\.([0-9]+).*", "\\1", aa_change_sim))]
sim_cds[, pos_in_codon := ((nt_pos_sim - 1) %% 3) + 1]

# confirm that this agrees with the snpEff uppercase annotation:
sim_cds[, codon_start_nt := ((nt_pos_sim - 1) %/% 3) * 3 + 1]
sim_cds[, codon_end_nt   := codon_start_nt + 2]

sim_cds[, snp_pos_in_codon := regexpr("[A-Z]", sub("/.*", "", codon_change_sim))]
sum(is.na(sim_cds$pos_in_codon) |
    is.na(sim_cds$snp_pos_in_codon) |
    sim_cds$pos_in_codon != sim_cds$snp_pos_in_codon
) # should be 0

# codon start/end in dsim3 coordinates
sim_cds[strand == "+",
    `:=`(
        dsim3_codon_start = pos_src - (pos_in_codon - 1),
        dsim3_codon_end   = pos_src + (3 - pos_in_codon)
    )]

sim_cds[strand == "-",
    `:=`(
        dsim3_codon_start = pos_src - (3 - pos_in_codon),
        dsim3_codon_end   = pos_src + (pos_in_codon - 1)
    )]

# codon start/end in dm6 coordinates
sim_cds[strand == "+", 
    `:=`(
        dm6_codon_start = pos - (pos_in_codon - 1),
        dm6_codon_end   = pos + (3 - pos_in_codon)
    )]

sim_cds[strand == "-", 
    `:=`(
        dm6_codon_start = pos - (3 - pos_in_codon),
        dm6_codon_end   = pos + (pos_in_codon - 1)
    )]

### check against reference
library(Biostrings)
library(data.table)

ref <- readDNAStringSet("/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna")
sim_cds[, codon_ref_annotation := toupper(sub("/.*", "", codon_change_sim))]
names(ref) <- trimws(names(ref))

sim_cds[, codon_genomic_raw := NA_character_]
for (chr in unique(sim_cds$sim_chr)) {
    idx <- which(sim_cds$sim_chr == chr)
    seq <- as.character(ref[[chr]])
    sim_cds$codon_genomic_raw[idx] <- substring(
        seq,
        sim_cds$dsim3_codon_start[idx],
        sim_cds$dsim3_codon_end[idx]
    )
}

sim_cds[strand == "+", codon_from_genome := codon_genomic_raw]
sim_cds[strand == "-", codon_from_genome := as.character( reverseComplement(DNAStringSet(as.character(codon_genomic_raw))))]

nrow(sim_cds[codon_from_genome == codon_ref_annotation])
nrow(sim_cds[codon_from_genome != codon_ref_annotation]) # should be 0. these I verified are at the CDS boundary 
sum(is.na(sim_cds$codon_from_genome)) # should be 0
sum(is.na(sim_cds$codon_ref_annotation)) # should be 0 

table(sim_cds$codon_from_genome == sim_cds$codon_ref_annotation)
View(sim_cds[codon_from_genome != codon_ref_annotation, .(cds_start, cds_end, dsim3_codon_start, dsim3_codon_end)]) # confirmed these anomalies overlap the CDS boundary
sim_cds <- sim_cds[codon_from_genome == codon_ref_annotation]
table(sim_cds$codon_from_genome == sim_cds$codon_ref_annotation)

# codon_genomic_start = genomic coord of 1st base on reference strand
#              _end = 3rd base 
# codon_genomic_raw = the actual codon (eg: CCT) read left to right on reference strand (+)
# codon_from_genome = genomic codon converted into 5' to 3' (reverse complmented for - strand)
# codon_ref_annotation = reference codon annotated by snpeff (from codon_change_sim)

### Dsim3 codon id: chr:[start]:[end]
sim_cds[, codon_id_dsim3 := paste(sim_chr, dsim3_codon_start, dsim3_codon_end, sep = ":" )]

### Dm6 codon id: chr:[start]:[end]
sim_cds[, codon_id_dm6_sim := paste(chr, dm6_codon_start, dm6_codon_end, sep = ":")]

# calculate number of snps per codon
sim_cds[, n_snps_in_codon_dsim3 := .N, by = codon_id_dsim3]
sim_cds[, n_snps_in_codon_dm6 := .N, by = codon_id_dm6_sim]

table(sim_cds$n_snps_in_codon_dsim3, useNA="ifany")
table(sim_cds$n_snps_in_codon_dm6, useNA="ifany")

sim_1snp_perCodon_dm6 <- copy(sim_cds[n_snps_in_codon_dm6 ==1])

sim_cds[, .N, by = .(chr, pos)][N > 1]
shared_dt <- merge(shared_dt, sim_cds[, .(chr, pos, codon_id_dsim3, codon_id_dm6_sim, strand, snp_pos_in_codon)], by=c("chr", "pos"), all.x=T)
setnames(shared_dt, c("strand", "snp_pos_in_codon"), c("strand_sim", "snp_pos_in_codon_sim"))

### 2B: melanogaster
# mel_dt <- readRDS("/scratch/ejy4bu/drosophila/DEST_remake/snpDT/mel/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.snp_dt_SynMissense.rds")
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

cds_gff <- gff[feature == "CDS" & !is.na(chr_dm6), .(chr_dm6, start, end, strand, phase, attributes)]

# extract genes and transcripts
cds_gff[, transcript_id := str_extract(attributes, "(?<=Parent=rna-)[^;,]+")]
mrna_map <- gff[
    feature == "mRNA",
    .(
        transcript_id = str_extract(attributes, "(?<=ID=rna-)[^;]+"),
        transcript_id_fbtr = str_extract(attributes, "FBtr[0-9]+")
    )
]
cds_gff[mrna_map, on = "transcript_id", transcript_id_fbtr := i.transcript_id_fbtr]

cds_gff[, gene_name := str_extract(attributes, "(?<=;gene=)[^;]+")]
cds_gff[, gene_id_fbgn := str_extract(attributes, "FBgn[0-9]+")]
cds_gff[, .N, by = .(chr_dm6, strand)]

# x prefix = cds_gff
# i prefix = shared_dt
mel_cds <- cds_gff[
    shared_dt[!is.na(variant.id_mel)], 
    on = .(
        chr_dm6 = chr,
        transcript_id_fbtr = transcript_id_mel,
        start <= pos, 
        end >= pos
    ), nomatch = 0,
    .(# GFF information
        variant.id_mel = i.variant.id_mel,
        chr = i.chr,
        pos = i.pos,
        ref_mel = i.ref_mel,
        alt_mel = i.alt_mel,

        effect_mel = i.effect_mel,
        codon_change_mel = i.codon_change_mel,
        aa_change_mel = i.aa_change_mel,
        aa_pos_mel = i.aa_pos_mel,

        variant.id_sim = i.variant.id_sim,
        ref_sim = i.ref_sim,
        alt_sim = i.alt_sim,
        effect_sim = i.effect_sim,
        codon_change_sim = i.codon_change_sim,
        aa_change_sim = i.aa_change_sim,
        aa_pos_sim = i.aa_pos_sim,

        ### gff info:
        chr_fromGFF = x.chr_dm6,
        cds_start = x.start,
        cds_end = x.end,
        strand = x.strand,
        phase = x.phase,
        attributes = x.attributes,
        transcript_id_fbtr = x.transcript_id_fbtr,
        transcript_id = x.transcript_id,
        gene_name = x.gene_name,
        gene_id_fbgn_fromGFF = x.gene_id_fbgn,
        gene_id_fbgn = i.gene_id_fbgn, 

        # nlp / allele frequency info
        global_af_mel = i.global_af_mel,
        poly_af_mel = i.poly_af_mel,
        nLocales_poly_mel = i.nLocales_poly_mel,
        nLocales_total_mel = i.nLocales_total_mel, 
        nLocales_missing_mel = i.nLocales_missing_mel,
        nLocales_poly_proj_round_mel = i.nLocales_poly_proj_round_mel,
        nLocales_poly_proj_binomial_mel = i.nLocales_poly_proj_binomial_mel

    )
]

mel_cds[, nt_pos_mel := as.integer(sub(".*c\\.([0-9]+).*", "\\1", aa_change_mel))]
mel_cds[, pos_in_codon := ((nt_pos_mel - 1) %% 3) + 1]
# confirm that this agrees with the snpEff uppercase annotation:
mel_cds[, codon_start_nt := ((nt_pos_mel - 1) %/% 3) * 3 + 1]
mel_cds[, codon_end_nt   := codon_start_nt + 2]

mel_cds[, snp_pos_in_codon := regexpr("[A-Z]", sub("/.*", "", codon_change_mel))]
sum(is.na(mel_cds$pos_in_codon) |
    is.na(mel_cds$snp_pos_in_codon) |
    mel_cds$pos_in_codon != mel_cds$snp_pos_in_codon
) # should be 0.

mel_cds[strand == "+",
    `:=`(
        codon_genomic_start = pos - (pos_in_codon - 1),
        codon_genomic_end   = pos + (3 - pos_in_codon)
    )]

mel_cds[strand == "-",
    `:=`(
        codon_genomic_start = pos - (3 - pos_in_codon),
        codon_genomic_end   = pos + (pos_in_codon - 1)
    )]

### check against reference
library(Biostrings)
library(data.table)

ref <- readDNAStringSet("/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna")
mel_cds[, codon_ref_annotation := toupper(sub("/.*", "", codon_change_mel))]
mel_cds[, codon_genomic_raw := NA_character_]

ref_chr <- ref[trimws(names(ref)) %in% chr_map$chr]

names(ref_chr) <- trimws(names(ref_chr))
for (chr in unique(mel_cds$chr)) {
    idx <- which(mel_cds$chr == chr)
    seq <- as.character(ref_chr[[chr]])

    mel_cds$codon_genomic_raw[idx] <- substring(
        seq,
        mel_cds$codon_genomic_start[idx],
        mel_cds$codon_genomic_end[idx]
    )
}


mel_cds[strand == "+", codon_from_genome := codon_genomic_raw]
mel_cds[strand == "-", codon_from_genome := as.character( reverseComplement(DNAStringSet(as.character(codon_genomic_raw))))]

nrow(mel_cds[codon_from_genome == codon_ref_annotation])
nrow(mel_cds[codon_from_genome != codon_ref_annotation]) # should be 0
sum(is.na(mel_cds$codon_from_genome)) # should be 0
sum(is.na(mel_cds$codon_ref_annotation)) # should be 0
table(mel_cds$codon_from_genome == mel_cds$codon_ref_annotation)
View(mel_cds[codon_from_genome != codon_ref_annotation, .(cds_start, cds_end, codon_genomic_start, codon_genomic_end)]) # confirmed these anomalies overlap the CDS boundary

# remove those codons that don't match bc they are at the CDS boundary
mel_cds <- mel_cds[codon_from_genome == codon_ref_annotation]
table(mel_cds$codon_from_genome == mel_cds$codon_ref_annotation)

# codon_genomic_start = genomic coord of 1st base on reference strand
#              _end = 3rd base 
# codon_genomic_raw = the actual codon (eg: CCT) read left to right on reference strand (+)
# codon_from_genome = genomic codon converted into 5' to 3' (reverse complmented for - strand)
# codon_ref_annotation = reference codon annotated by snpeff (from codon_change_sim)


### dm6 codon id: chr:[start]:[end]
mel_cds[, codon_id_dm6_mel := paste(chr, codon_genomic_start, codon_genomic_end, sep = ":")]

mel_cds[, n_snps_in_codon := .N, by = codon_id_dm6_mel]
mel_1snp_perCodon <- copy(mel_cds[n_snps_in_codon ==1])


mel_cds[, .N, by = .(chr, pos)][N > 1]
shared_dt <- merge(shared_dt, mel_cds[, .(chr, pos, gene_id_fbgn_fromGFF, codon_id_dm6_mel, strand, snp_pos_in_codon)], by=c("chr", "pos"), all.x=T)
setnames(shared_dt, c("strand", "snp_pos_in_codon"), c("strand_mel", "snp_pos_in_codon_mel"))

###############################################################################################
### 3A: D. simulans: Annotate shared table for 1 snp per codon
###############################################################################################
shared_dt[, keep_sim := fifelse(
    is.na(variant.id_sim), "No_sim_variant",
    NA_character_
)]
table(shared_dt$keep_sim, useNA="ifany")

shared_dt[!is.na(variant.id_sim), keep_sim := fifelse(
    variant.id_sim %in% sim_1snp_perCodon_dm6$variant.id_sim, 
    "TRUE", 
    "Multiple_SNPs_perRefCodon_dm6"
)]
table(shared_dt$keep_sim, useNA="ifany")


###############################################################################################
### 3B: D. melanogaster: Annotated shared table for 1 snp per codon
###############################################################################################
shared_dt[, keep_mel := fifelse(
    is.na(variant.id_mel), "No_mel_variant",
    NA_character_
)]
table(shared_dt$keep_mel, useNA="ifany")

shared_dt[!is.na(variant.id_mel), keep_mel := fifelse(
    variant.id_mel %in% mel_1snp_perCodon$variant.id_mel, 
    "TRUE", 
    "Multiple_SNPs_perRefCodon"
)]
table(shared_dt$keep_mel, useNA="ifany")


###############################################################################################
### 4: Filter for unassigned transcripts / codons missing from CDS 
###############################################################################################

### 4A: simulans
missing_from_sim_cds <- shared_dt[!is.na(variant.id_sim) & !(variant.id_sim %in% sim_cds$variant.id_sim)]
nrow(missing_from_sim_cds) # they are all unassigned transcripts.

shared_dt[variant.id_sim %in% missing_from_sim_cds$variant.id_sim, keep_sim := "Missing_fromCDS"]
table(shared_dt$keep_sim, useNA="ifany")

### 4B: melanogaster 
missing_from_mel_cds <- shared_dt[!is.na(variant.id_mel) & !(variant.id_mel %in% mel_cds$variant.id_mel)]
nrow(missing_from_mel_cds)

shared_dt[variant.id_mel %in% missing_from_mel_cds$variant.id_mel, keep_mel := "Missing_fromCDS"]
table(shared_dt$keep_mel, useNA="ifany")


###############################################################################################
### 5: Flag codons mapped 2 to 1
###############################################################################################

dup_dm6_codons <- sim_1snp_perCodon_dm6[, .N, by = codon_id_dm6_sim][N > 1]

shared_dt[variant.id_sim %in% sim_1snp_perCodon_dm6[
        codon_id_dm6_sim %in% dup_dm6_codons$codon_id_dm6_sim,
        variant.id_sim
    ],
    keep_sim := "MAP2TO1"
]
# none! already removed when building sim_1snp_perCodon_dm6 
    # because 2 dsim3 codons mapping to 1 dm6 codon would be counted as 2 snps in that codon
table(shared_dt$keep_sim, useNA="ifany")

###############################################################################################
### 6: Identify codons with shifted reading frame
###############################################################################################

mel_codons <- unique(
    shared_dt[
        !is.na(variant.id_mel),
        .(
            codon_id = codon_id_dm6_mel,
            chr = chr,
            start = as.integer(tstrsplit(codon_id_dm6_mel, ":", fixed = TRUE)[[2]]),
            end = as.integer(tstrsplit(codon_id_dm6_mel, ":", fixed = TRUE)[[3]])
        )
    ]
)
sim_codons <- unique(
    shared_dt[
        !is.na(variant.id_sim),
        .(
            codon_id = codon_id_dm6_sim,
            chr = chr,
            start = as.integer(tstrsplit(codon_id_dm6_sim, ":", fixed = TRUE)[[2]]),
            end = as.integer(tstrsplit(codon_id_dm6_sim, ":", fixed = TRUE)[[3]])
        )
    ]
)

overlap <- mel_codons[sim_codons, on = .(chr, start <= end, end >= start),
    nomatch = 0,
    .(codon_id_mel=x.codon_id, codon_id_sim=i.codon_id),
    allow.cartesian = TRUE
]

shifted_overlap <- overlap[codon_id_mel != codon_id_sim]

shared_dt[
    codon_id_dm6_mel %in% shifted_overlap$codon_id_mel,
    keep_mel := "shifted_readingFrame"
]
table(shared_dt$keep_mel, useNA="ifany")

shared_dt[
    codon_id_dm6_sim %in% shifted_overlap$codon_id_sim,
    keep_sim := "shifted_readingFrame"
]
table(shared_dt$keep_sim, useNA="ifany")

###############################################################################################
### 7: Identify codons on opposite strands
###############################################################################################

# nrow(shared_dt[!is.na(variant.id_mel) & !is.na(variant.id_sim) & 
#     codon_id_dm6_mel != codon_id_dm6_sim & 
#     keep_sim == "TRUE" & keep_mel == "TRUE"]) # should be 0 
shared_dt[, codon_id_dm6 := fifelse(
    (!is.na(variant.id_mel) & keep_mel=="TRUE"), codon_id_dm6_mel, 
    fifelse(
        (!is.na(variant.id_sim) & keep_sim=="TRUE"), codon_id_dm6_sim, NA_character_))]

# !is.na(codon_id_dm6_mel) & !is.na(codon_id_dm6_sim)
# opp_codons <- shared_dt[,
#     .(
#         has_plus  = any(strand_mel == "+" | strand_sim == "+", na.rm = TRUE),
#         has_minus = any(strand_mel == "-" | strand_sim == "-", na.rm = TRUE)
#     ),
#     by = codon_id_dm6
# ][has_plus & has_minus, codon_id_dm6]

opp_codons <- shared_dt[
    !is.na(variant.id_mel) & !is.na(variant.id_sim) &
    # keep_mel == "TRUE" & keep_sim == "TRUE" &
    # !is.na(strand_mel) & !is.na(strand_sim) &
    strand_mel != strand_sim,
    unique(codon_id_dm6)
]

shared_dt[!is.na(codon_id_dm6) & codon_id_dm6%in%opp_codons, `:=` (keep_mel = "oppStrand", keep_sim ="oppStrand")]
table(shared_dt$keep_sim, useNA="ifany")
table(shared_dt$keep_mel, useNA="ifany")

saveRDS(shared_dt, paste0(out_dir, "dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.annotatedSim.annotatedMel.rds"))


###############################################################################################
### 8: Clean up dt
###############################################################################################
filtered_dt <- shared_dt[!is.na(variant.id_mel) | !is.na(variant.id_sim)]

# filtered_dt[, maf := NULL]

setnames(filtered_dt, c("chr_src", "pos_src", "ref_sim", "alt_sim", "ref_src", "alt_src"), 
    c("chr_dsim3", "pos_dsim3", "ref_sim_dm6", "alt_sim_dm6", "ref_sim_dsim3", "alt_sim_dsim3"))
mel_cols <- c(
    "ref_mel", "alt_mel", 
    "codon_id_dm6_mel", 
    "codon_ref_mel", "codon_alt_mel",
    "snp_pos_in_codon_mel", 
    "aa_pos_mel", 
    "aa_ref_mel", "aa_alt_mel", 
    "gene_mel", "gene_id_fbgn_fromGFF", 
    "transcript_id_mel", 
    "effect_mel", 
    "strand_mel", 
    "global_af_mel", "poly_af_mel", "poly_maf_mel",
    "nLocales_poly_mel", "nLocales_total_mel", "nLocales_missing_mel", 
    "nLocales_poly_proj_round_mel", "nLocales_poly_proj_binomial_mel"

)
### sim cols:
sim_cols <- c(
    "chr_dsim3", "pos_dsim3", 
    "ref_sim_dm6", "alt_sim_dm6", "ref_sim_dsim3", "alt_sim_dsim3",
    "codon_id_dm6_sim", "codon_id_dsim3",
    "codon_ref_sim", "codon_alt_sim",
    "snp_pos_in_codon_sim", 
    "aa_pos_sim",
    "aa_ref_sim", "aa_alt_sim",
    "gene_sim",
    "transcript_id_sim",
    "effect_sim", 
    "flip", "swap", 
    "strand_sim", 
    "global_af_sim", "poly_af_sim", "poly_maf_sim",
    "nLocales_poly_sim", "nLocales_total_sim", "nLocales_missing_sim", 
    "nLocales_poly_proj_round_sim", "nLocales_poly_proj_binomial_sim"
)

filtered_dt[
    is.na(variant.id_mel),
    (mel_cols) := lapply(.SD, function(x) {
        if (is.character(x)) NA_character_ else NA
    }),
    .SDcols = mel_cols
]

filtered_dt[
    is.na(variant.id_sim),
    (sim_cols) := lapply(.SD, function(x) {
        if (is.character(x)) NA_character_ else NA
    }),
    .SDcols = sim_cols
]

new_order <- c(
    # meta
    "variant.id_mel", "variant.id_sim", 
    "chr", "pos", 
    "chr_dsim3", "pos_dsim3", 
    # nucleotides
    "ref_mel", "alt_mel",
    "ref_sim_dm6", "alt_sim_dm6",
    "ref_sim_dsim3", "alt_sim_dsim3",
    "snp_pos_in_codon_mel", "snp_pos_in_codon_sim",
    # codons
    "codon_id_dm6", 
    "codon_id_dm6_mel", "codon_id_dm6_sim", 
    "codon_id_dsim3",
    "codon_ref_mel", "codon_alt_mel",
    "codon_ref_sim", "codon_alt_sim",
    # amino acids
    "aa_pos_mel", "aa_ref_mel", "aa_alt_mel", 
    "aa_pos_sim", "aa_ref_sim", "aa_alt_sim", 
    # gene info
    "gene_mel", "gene_sim", "gene_id_fbgn", "gene_id_fbgn_fromGFF",
    "transcript_id_mel", "transcript_id_sim", 
    # nlp / AF info 
    "global_af_mel", "poly_af_mel", 
    "global_af_sim", "poly_af_sim",
    
    "nLocales_poly_mel", "nLocales_total_mel", "nLocales_missing_mel", 
    "nLocales_poly_proj_round_mel", "nLocales_poly_proj_binomial_mel", 
    "nLocales_poly_sim", "nLocales_total_sim", "nLocales_missing_sim", 
    "nLocales_poly_proj_round_sim", "nLocales_poly_proj_binomial_sim",

    "effect_mel", "effect_sim",
    "flip", "swap", "strand_sim", "strand_mel", 
    "keep_mel", "keep_sim"
)

new_order <- intersect(new_order, names(filtered_dt))
filtered_dt <- filtered_dt[, ..new_order]
saveRDS(filtered_dt, paste0(out_dir, "dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.annotatedSim.annotatedMel.cleaner.rds"))

### end new script 








# # shared_dt[, keep_sim := as.character(variant.id_sim %in% sim_1snp_perCodon$variant_id_sim)]

# # shared_dt[, keep_sim := fifelse(
# #     is.na(variant.id_sim),
# #     NA_character_,
# #     as.character(variant.id_sim %in% sim_1snp_perCodon$variant_id_sim)
# # )]

# # this annotates all rows TRUE or FALSE depending on if they are in the sim dt 
#     # so FALSE includes rows that didn't pass multiple snp per codon filter AND 
#     #   those that have no sim annotation at that row
# # the number of TRUE rows == number of rows in sim_1snp_perCodon



# # missing_from_sim_cds <- shared_dt[
# #     keep_sim == FALSE & !is.na(variant.id_sim) &
# #     !(variant.id_sim %in% sim_cds$variant_id_sim)
# # ]

# # nrow(missing_from_sim_cds)
# # they are all unassigned transcripts.. 
# # except for 1, which i removed because it crossed the CDS boundary so snpeff codon an ref fasta codon didn't match...

# shared_dt[variant.id_sim %in% missing_from_sim_cds$variant.id_sim, keep_sim := "UNASSIGNED"]
# table(shared_dt$keep_sim, useNA="ifany")

# # saveRDS(shared_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.rds")
# # saveRDS(shared_dt, paste0(out_dir, "dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.merge_unfilt.annotatedSim.rds"))
# saveRDS(shared_dt, "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/shared/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.annotatedSim.rds")

# ############################################################################################
# ### repeat for mel: 1 snp per codon



# # there are 6197 rows that are missing from mel cds (annotated as FALSE but have variant id mel)
# missing_from_mel_cds <- shared_dt[
#     keep_mel == FALSE & !is.na(variant.id_mel) &
#     !(variant.id_mel %in% mel_cds$variant_id_mel)
# ]

# nrow(missing_from_mel_cds)
# nrow(missing_from_mel_cds[transcript_id_mel%like%"unassigned"]) # 0? 
# shared_dt[
#     !is.na(variant.id_mel) &
#     !(variant.id_mel %in% mel_cds$variant_id_mel),
#     keep_mel := "MISSING_CDS"
# ]

# # shared_dt[variant.id_mel %in% missing_from_mel_cds$variant.id_mel, keep_mel := "UNASSIGNED"]
# table(shared_dt$keep_mel, useNA="ifany")

# # saveRDS(shared_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.rds")
# # saveRDS(shared_dt, paste0(out_dir, "dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.merge_unfilt.annotatedSim.annotatedMel.rds"))
# saveRDS(shared_dt, "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/shared/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.annotatedSim.annotatedMel.rds")


# ############################################################################################
# ### check if any of the sim liftover results show up as in the wrong spot for mel 
#     # in order to determine same / different codon snps

# ### add dm6 codon id for sim variants:
# sim_1snp_perCodon[
#     strand == "+",
#     `:=`(
#         dm6_codon_start = pos - (snp_pos_in_codon - 1),
#         dm6_codon_end   = pos + (3 - snp_pos_in_codon)
#     )
# ]

# sim_1snp_perCodon[
#     strand == "-",
#     `:=`(
#         dm6_codon_start = pos - (3 - snp_pos_in_codon),
#         dm6_codon_end   = pos + (snp_pos_in_codon - 1)
#     )
# ]

# sim_1snp_perCodon[
#     ,
#     dm6_codon_id := paste(
#         chr,
#         dm6_codon_start,
#         dm6_codon_end,
#         sep = ":"
#     )
# ]

# ### rename mel codon id for dm6 coordinates
# mel_1snp_perCodon[
#     ,
#     dm6_codon_id := codon_id
# ]

# ### ID codons where 2 sim dsim3 codons map to 1 dm6 codon
# dup_dm6_codons <- sim_1snp_perCodon[
#     ,
#     .N,
#     by = dm6_codon_id
# ][N > 1]


# shared_dt[
#     variant.id_sim %in% sim_1snp_perCodon[
#         dm6_codon_id %in% dup_dm6_codons$dm6_codon_id,
#         variant_id_sim
#     ],
#     keep_sim := "MAP2TO1"
# ]

# table(shared_dt$keep_sim, useNA = "ifany")
# shared_dt_all <- shared_dt[
#     sim_1snp_perCodon,
#     on = .(variant.id_sim = variant_id_sim),
#     `:=`(
#         strand_sim = i.strand,
#         snp_pos_in_codon_sim = i.snp_pos_in_codon,
#         codon_id_dsim3 = i.codon_id,
#         codon_id_dm6_sim = i.dm6_codon_id
#     )
# ]

# shared_dt_all <- shared_dt_all[
#     mel_1snp_perCodon,
#     on = .(variant.id_mel = variant_id_mel),
#     `:=`(
#         strand_mel = i.strand,
#         snp_pos_in_codon_mel = i.snp_pos_in_codon,
#         codon_id_dm6_mel = i.dm6_codon_id,
#         gene_id_fbgn_fromGFF=i.gene_id_fbgn_fromGFF, 
#         gene_id_fbgn=i.gene_id_fbgn
#     )
# ]
    
# # (, .(strand, snp_pos_in_codon, codon_id, dm6_codon_id),)]
# shared_dt_all[, classification := NA_character_] # add empty classification colun

# shared_dt_cleaner <- shared_dt_all[, .(
#     ### meta stuff 
#     classification,  
#     variant.id_mel, variant.id_sim,
#     chr_dm6=chr, pos_dm6=pos,
#     chr_dsim3=chr_src, pos_dsim3=pos_src,
    
#     ### nt stuff:
#     ref_mel, alt_mel, 
#     ref_sim_dm6=ref_sim, alt_sim_dm6=alt_sim, 
#     ref_sim_dsim3=ref_src, alt_sim_dsim3=alt_src, 

#     ### codon stuff:
#     codon_id_dm6_mel, codon_id_dm6_sim, codon_id_dsim3,
#     codon_ref_mel, codon_alt_mel, codon_ref_sim, codon_alt_sim,
#     snp_pos_in_codon_mel, snp_pos_in_codon_sim,

#     ### amino acid stuff:
#     aa_pos_mel, aa_pos_sim,
#     aa_ref_mel, aa_alt_mel, aa_ref_sim, aa_alt_sim,

#     ### gene / transcript stuff:
#     gene_mel, gene_sim, 
#     gene_id_fbgn,
#     gene_id_fbgn_fromGFF,

#     transcript_id_mel, transcript_id_sim,

#     ### frequency stuff:
#     # n_samps_mel, n_samps_sim, 
#     # af_mel, maf_mel, af_sim, maf_sim, 
#     effect_mel, effect_sim, 

#     ### extras:
#     flip, swap, 
#     strand_mel, strand_sim,
#     keep_mel, keep_sim
# )]

# # saveRDS(shared_dt_cleaner, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.cleaner.rds")
# saveRDS(shared_dt_cleaner, "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/shared/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.annotatedSim.annotatedMel.cleaner.rds")
# ### clean up rows where keep_mel / sim != TRUE

# filtered_dt <- shared_dt_cleaner
# filtered_dt[keep_mel != "TRUE", variant.id_mel := NA_integer_]
# filtered_dt[keep_sim != "TRUE", variant.id_sim := NA_integer_]

# filtered_dt <- filtered_dt[!is.na(variant.id_mel) | !is.na(variant.id_sim)]
# ### mel cols:
# mel_cols <- c(
#     "ref_mel", "alt_mel", 
#     "codon_id_dm6_mel", 
#     "codon_ref_mel", "codon_alt_mel",
#     "snp_pos_in_codon_mel", 
#     "aa_pos_mel", 
#     "aa_ref_mel", "aa_alt_mel", 
#     "gene_mel", "gene_id_fbgn_fromGFF", 
#     "transcript_id_mel", 
#     # "n_samps_mel", 
#     # "af_mel", "maf_mel", 
#     "effect_mel", 
#     "strand_mel"
# )
# ### sim cols:
# sim_cols <- c(
#     "chr_dsim3", "pos_dsim3", 
#     "ref_sim_dm6", "alt_sim_dm6", "ref_sim_dsim3", "alt_sim_dsim3",
#     "codon_id_dm6_sim", "codon_id_dsim3",
#     "codon_ref_sim", "codon_alt_sim",
#     "snp_pos_in_codon_sim", 
#     "aa_pos_sim",
#     "aa_ref_sim", "aa_alt_sim",
#     "gene_sim",
#     "transcript_id_sim",
#     # "n_samps_sim", 
#     # "af_sim", "maf_sim",
#     "effect_sim", 
#     "flip", "swap", 
#     "strand_sim"
# )

# filtered_dt[
#     is.na(variant.id_mel),
#     (mel_cols) := lapply(.SD, function(x) {
#         if (is.character(x)) NA_character_ else NA
#     }),
#     .SDcols = mel_cols
# ]

# filtered_dt[
#     is.na(variant.id_sim),
#     (sim_cols) := lapply(.SD, function(x) {
#         if (is.character(x)) NA_character_ else NA
#     }),
#     .SDcols = sim_cols
# ]

# # saveRDS(filtered_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.cleaner.1snpCodon.rds")


# ### working off filtered_dt


# ### checking if there are codons that overlap (different reading frame)
# mel_codons <- unique(
#     filtered_dt[
#         !is.na(variant.id_mel),
#         .(
#             codon_id = codon_id_dm6_mel,
#             chr = chr_dm6,
#             start = as.integer(tstrsplit(codon_id_dm6_mel, ":", fixed = TRUE)[[2]]),
#             end = as.integer(tstrsplit(codon_id_dm6_mel, ":", fixed = TRUE)[[3]])
#         )
#     ]
# )

# sim_codons <- unique(
#     filtered_dt[
#         !is.na(variant.id_sim),
#         .(
#             codon_id = codon_id_dm6_sim,
#             chr = chr_dm6,
#             start = as.integer(tstrsplit(codon_id_dm6_sim, ":", fixed = TRUE)[[2]]),
#             end = as.integer(tstrsplit(codon_id_dm6_sim, ":", fixed = TRUE)[[3]])
#         )
#     ]
# )

# overlap <- mel_codons[
#     sim_codons,
#     on = .(
#         chr,
#         start <= end,
#         end >= start
#     ),
#     nomatch = 0,
#     .(
#         codon_id_mel=x.codon_id,
#         codon_id_sim=i.codon_id
#     ),
#     allow.cartesian = TRUE
# ]

# shifted_overlap <- overlap[
#     codon_id_mel != codon_id_sim
# ]


# filtered_dt[
#     codon_id_dm6_mel %in% shifted_overlap$codon_id_mel,
#     keep_mel := "OVERLAP"
# ]

# filtered_dt[
#     codon_id_dm6_sim %in% shifted_overlap$codon_id_sim,
#     keep_sim := "OVERLAP"
# ]



# ### inspecting # snps per codons

# mel_pos <- filtered_dt[
#     !is.na(variant.id_mel),
#     .(codon_id = codon_id_dm6_mel,
#       mel_pos = snp_pos_in_codon_mel)
# ]

# sim_pos <- filtered_dt[
#     !is.na(variant.id_sim),
#     .(codon_id = codon_id_dm6_sim,
#       sim_pos = snp_pos_in_codon_sim)
# ]

# shared <- merge(mel_pos, sim_pos, by = "codon_id")

# shared[, same_position := mel_pos == sim_pos]

# table(shared$same_position,useNA="ifany")

# shared[, tsp_samePos := same_position]

# tsp_map <- unique(shared[, .(
#     codon_id = codon_id,
#     tsp_samePos
# )])

# filtered_dt[
#     tsp_map,
#     on = .(
#         codon_id_dm6_mel = codon_id
#     ),
#     tsp_samePos := i.tsp_samePos
# ]

# filtered_dt[
#     tsp_map,
#     on = .(
#         codon_id_dm6_sim = codon_id
#     ),
#     tsp_samePos := i.tsp_samePos
# ]

# table(filtered_dt$tsp_samePos, useNA="ifany") 
# # TRUE = number of variants with same position variant
# # FALSE = 2x number of codons with diff position snps
# # NA = species specific variants OR codon id did not match (reading frame was off)


# # saveRDS(filtered_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.cleaner.1snpCodon.rds")
# # saveRDS(filtered_dt, paste0(out_dir, "dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.merge_unfilt.annotatedSim.annotatedMel.cleaner.1snpCodon.rds"))
# saveRDS(filtered_dt, "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/shared/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.annotatedSim.annotatedMel.cleaner.1snpCodon.rds")






# ### from pre-classification script - to integrate
# filtered_dt[, codon_id_dm6  := fifelse(
#     !is.na(codon_id_dm6_mel),
#     codon_id_dm6_mel,
#     codon_id_dm6_sim)]

# # realized I filtered for opp codons earlier but missed some... 
# opp_codons <- filtered_dt[
#     ,
#     .(
#         has_plus  = any(strand_mel == "+" | strand_sim == "+", na.rm = TRUE),
#         has_minus = any(strand_mel == "-" | strand_sim == "-", na.rm = TRUE)
#     ),
#     by = codon_id_dm6
# ][has_plus & has_minus, codon_id_dm6]

# filtered_dt[codon_id_dm6%in%opp_codons, `:=` (keep_mel = "oppStrnd", keep_sim ="oppStrnd")]

# # opp_strand <- filtered_dt[
# #     keep_mel=="TRUE" & keep_sim=="TRUE" &
# #     !is.na(codon_id_dm6_mel) &
# #     !is.na(codon_id_dm6_sim) &
# #     !is.na(strand_mel) &
# #     !is.na(strand_sim) &
# #     strand_mel != strand_sim
# # ]
# # filtered_dt[codon_id_dm6_mel%in%opp_strand$codon_id_dm6_mel, `:=` (keep_mel="oppStrnd", keep_sim="oppStrnd")]


# filtered_dt[keep_mel=="OVERLAP", samePos := "OVERLAP"]
# filtered_dt[keep_sim=="OVERLAP", samePos := "OVERLAP"]


# ### mel cols:
# mel_cols <- c(
#     "variant.id_mel",
#     "ref_mel", "alt_mel", 
#     "codon_id_dm6_mel", 
#     "codon_ref_mel", "codon_alt_mel",
#     "snp_pos_in_codon_mel", 
#     "aa_pos_mel", 
#     "aa_ref_mel", "aa_alt_mel", 
#     "gene_mel", "gene_id_fbgn_fromGFF", 
#     "transcript_id_mel", 
#     "n_samps_mel", 
#     "af_mel", #"maf_mel", 
#     "effect_mel", 
#     "strand_mel"
# )
# ### sim cols:
# sim_cols <- c(
#     "variant.id_sim",
#     "chr_dsim3", "pos_dsim3", 
#     "ref_sim_dm6", "alt_sim_dm6", "ref_sim_dsim3", "alt_sim_dsim3",
#     "codon_id_dm6_sim", "codon_id_dsim3",
#     "codon_ref_sim", "codon_alt_sim",
#     "snp_pos_in_codon_sim", 
#     "aa_pos_sim",
#     "aa_ref_sim", "aa_alt_sim",
#     "gene_sim",
#     "transcript_id_sim",
#     "n_samps_sim", 
#     "af_sim", # "maf_sim",
#     "effect_sim", 
#     "flip", "swap", 
#     "strand_sim"
# )
# filtered_dt[keep_mel != "TRUE", (mel_cols) := NA]
# filtered_dt[keep_sim != "TRUE", (sim_cols) := NA]

# # filtered_dt[keep_mel != "TRUE", variant.id_mel := NA_character_]
# new_dt <- filtered_dt[!is.na(variant.id_mel) | !is.na(variant.id_sim)]
# # new_dt <- new_dt[keep_mel == "TRUE" | is.na(keep_mel)]
# # new_dt <- new_dt[keep_sim == "TRUE" | is.na(keep_sim)] 
# ## ^this was wrong because if mel variant was present and failed, the whole row would be removed before sim could be assessed
# # new_dt <- new_dt[keep_mel == "TRUE" | keep_sim == "TRUE"]


# new_dt <- new_dt[
#     (keep_mel == "TRUE" | is.na(keep_mel)) | 
#     (keep_sim == "TRUE" | is.na(keep_sim))
# ]

# new_dt[, tsp_samePos := snp_pos_in_codon_mel == snp_pos_in_codon_sim]



# new_dt[, .N, by = codon_id_dm6][, table(N)]
# new_dt[is.na(tsp_samePos), .N, by = codon_id_dm6][, table(N)] # mixture of diffPos snps and species-specific sites 
# new_dt[!is.na(tsp_samePos), .N, by = codon_id_dm6][, table(N)] # all the samePos sites (all TRUE), so all 1 codon only... good!

# # View(new_dt[is.na(tsp_samePos), .N, by = codon_id_dm6])
# # # get the list of these codons where only 1 codon id dm6 but is not same pos and remove... 

# # remove_codons <- new_dt[
# #     is.na(tsp_samePos),
# #     .N,
# #     by = codon_id_dm6
# # ][N == 1, codon_id_dm6]

# # length(remove_codons)
# # new_dt <- new_dt[!(codon_id_dm6 %in% remove_codons)]
# ### then classify!!!

# tsp_map <- new_dt[
#     ,
#     .(
#         pos_mel = snp_pos_in_codon_mel[!is.na(snp_pos_in_codon_mel)][1],
#         pos_sim = snp_pos_in_codon_sim[!is.na(snp_pos_in_codon_sim)][1]
#     ),
#     by = codon_id_dm6
# ]

# tsp_map[, samePos := pos_mel == pos_sim]

# new_dt[
#     tsp_map,
#     on = "codon_id_dm6",
#     tsp_samePos := i.samePos
# ]

# # check everything:
# table(tsp_map$samePos, useNA = "ifany")
# table(new_dt$tsp_samePos, useNA="ifany")


# new_dt[
#      ,
#      .(
#          N_rows = .N,
#          n_mel = sum(!is.na(variant.id_mel)),
#          n_sim = sum(!is.na(variant.id_sim))
#      ),
#      by = .(codon_id_dm6, tsp_samePos)
#  ][
#      ,
#      .N,
#      by = .(tsp_samePos, N_rows, n_mel, n_sim)
#  ][order(tsp_samePos, N_rows, n_mel, n_sim)]