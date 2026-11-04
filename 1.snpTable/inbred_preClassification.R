library(data.table)
library(stringr)
library(Biostrings)
library(SeqArray)

#######################################
### sim signor
# lifted_vcf <- "/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.dm6.sorted.m2.M2.vcf.gz"

# genome_fna <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"
# cds_fasta <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_cds_genomic.fna"

# sim_dt <- readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.snp_dt_SynMissense.rds")
# gds_file <- "/scratch/ejy4bu/drosophila/inbred/gds_files/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.dm6.sorted.m2.M2.goodSamps.goodSites.gds"

shared_dt <-readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.10_2026.rds")

# out_dir <- "/scratch/ejy4bu/drosophila/inbred/snpDT/simMap/"
# out_rds <- paste0(out_dir,  "dsim3.signor.all_quality_variants_codonMap.10_2026.rds")

###############################################################################################
### 1. recalculate AF:
###############################################################################################

# 1A: Recalculate AF for D. simulans: 
sim_gds <- seqOpen("/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.dm6.sorted.m2.M2.goodSamps.goodSites.gds")

snp_dt <- copy(shared_dt)
shared_dt[, `:=` (af_sim_old = af_sim, n_samps_sim_old = n_samps_sim, af_sim = NULL, n_samps_sim = NULL)]
snp_dt[, `:=`(af_sim_old  = af_sim, n_samps_sim_old = n_samps_sim)]
sim_candidates <- snp_dt[!is.na(variant.id_sim), .(chr, pos)]

seqResetFilter(sim_gds)
seqSetFilterPos(sim_gds, chr=sim_candidates$chr, pos=sim_candidates$pos)
sim_freq_dt <- data.table(
    chr = seqGetData(sim_gds, "chromosome"), 
    pos = seqGetData(sim_gds, "position"), 
    af_sim = seqAlleleFreq(sim_gds, ref.allele = 1L)
)

sim_freq_dt[, maf_sim := pmin(af_sim, 1-af_sim)]
snp_dt[sim_freq_dt, on = .(chr, pos), `:=` (af_sim = i.af_sim, maf_sim = i.maf_sim)]
snp_dt[is.na(variant.id_sim), `:=` (af_sim = NA_real_, maf_sim = NA_real_)]

shared_dt <- merge(shared_dt, snp_dt[, .(chr, pos, af_sim, maf_sim, n_samps_sim)], by=c("chr", "pos"), all.x=T)
nrow(shared_dt[!is.na(variant.id_sim) & af_sim > af_sim_old])

seqClose(sim_gds)

# 1B: Recalculate AF for D. melanogaster:
mel_gds <- seqOpen("/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.goodSamps.goodSites.gds")

snp_dt <- copy(shared_dt)
shared_dt[, `:=` (af_mel_old = af_mel, n_samps_mel_old = n_samps_mel, af_mel = NULL, n_samps_mel = NULL)]
snp_dt[, `:=`(af_mel_old  = af_mel, n_samps_mel_old = n_samps_mel)]
mel_candidates <- snp_dt[!is.na(variant.id_mel), .(chr, pos)]

seqResetFilter(mel_gds)
seqSetFilterPos(mel_gds, chr=mel_candidates$chr, pos=mel_candidates$pos)
mel_freq_dt <- data.table(
    chr = seqGetData(mel_gds, "chromosome"), 
    pos = seqGetData(mel_gds, "position"), 
    af_mel = seqAlleleFreq(mel_gds, ref.allele = 1L)
)
mel_freq_dt[, maf_mel := pmin(af_mel, 1-af_mel)]
snp_dt[mel_freq_dt, on = .(chr, pos), `:=` (af_mel = i.af_mel, maf_mel = i.maf_mel)]
snp_dt[is.na(variant.id_mel), `:=` (af_mel = NA_real_, maf_mel = NA_real_)]
shared_dt <- merge(shared_dt, snp_dt[, .(chr, pos, af_mel, maf_mel, n_samps_mel)], by=c("chr", "pos"), all.x=T)
nrow(shared_dt[!is.na(variant.id_mel) & af_mel > af_mel_old])
seqClose(mel_gds)

# safe copy:
shared_dt_noMAFfilter <- copy(shared_dt)
nrow(shared_dt_noMAFfilter)

saveRDS(shared_dt_noMAFfilter, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.newAF.noMAFfilter.rds")

shared_dt[!is.na(variant.id_mel) & (maf_mel < 0.05 | maf_mel > 0.95),
    variant.id_mel := NA
]
shared_dt[
    !is.na(variant.id_sim) & (maf_sim < 0.05 | maf_sim > 0.95),
    variant.id_sim := NA
]

shared_dt <- shared_dt[!is.na(variant.id_sim) | !is.na(variant.id_mel), ]
nrow(shared_dt)

saveRDS(shared_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.newAF.0.05MAFfilter.rds")

###############################################################################################
### 2. Extract GFF info: CDS regions; chr, start, end, strand, phase 
###############################################################################################
### 2A: simulans
sim_gff <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.gff"

## gff format:
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
    shared_dt[!is.na(variant.id_sim) & af_sim >= 0.05 & af_sim <= 0.95], 
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
        af_mel = i.af_mel,
        effect_mel = i.effect_mel,
        codon_change_mel = i.codon_change_mel,
        aa_change_mel = i.aa_change_mel,
        aa_pos_mel = i.aa_pos_mel,

        variant.id_sim = i.variant.id_sim,
        ref_sim = i.ref_sim,
        alt_sim = i.alt_sim,
        af_sim = i.af_sim,
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
        swap = i.swap
    )
]

sim_cds[, nt_pos_sim := as.integer(sub(".*c\\.([0-9]+).*", "\\1", aa_change_sim))]
sim_cds[, pos_in_codon := ((nt_pos_sim - 1) %% 3) + 1]

# confirm that this agrees with the snpEff annotation:
sim_cds[, codon_start_nt := ((nt_pos_sim - 1) %/% 3) * 3 + 1]
sim_cds[, codon_end_nt   := codon_start_nt + 2]
sim_cds[, snp_pos_in_codon := regexpr("[A-Z]", sub("/.*", "", codon_change_sim))]
sum(is.na(sim_cds$pos_in_codon) |
    is.na(sim_cds$snp_pos_in_codon) |
    sim_cds$pos_in_codon != sim_cds$snp_pos_in_codon
) # should be 0

### assign strandedness 
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

### check against reference fasta 
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
### assign codons by strandedness 
sim_cds[strand == "+", codon_from_genome := codon_genomic_raw]
sim_cds[strand == "-", codon_from_genome := as.character(reverseComplement(DNAStringSet(as.character(codon_genomic_raw))))]

nrow(sim_cds[codon_from_genome == codon_ref_annotation])
nrow(sim_cds[codon_from_genome != codon_ref_annotation]) # should be 0
sum(is.na(sim_cds$codon_from_genome)) # should be 0
sum(is.na(sim_cds$codon_ref_annotation)) # should be 0 

table(sim_cds$codon_from_genome == sim_cds$codon_ref_annotation) 

# codon_from_genome matches codon_ref_annotation -- good!

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

shared_dt <- merge(shared_dt, sim_cds[, .(chr, pos, codon_id_dsim3, codon_id_dm6_sim, strand, snp_pos_in_codon)], by=c("chr", "pos"), all.x=T)
setnames(shared_dt, c("strand", "snp_pos_in_codon"), c("strand_sim", "snp_pos_in_codon_sim"))

### 2B: melanogaster
# mel_dt <- readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/DGRP2.source_BCM-HGSC.dm6.snp_dt_SynMissense.10_2026.rds")
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

# cds_gff[, transcript_id_fbtr := str_extract(attributes, "FBtr[0-9]+")]
cds_gff[, gene_name := str_extract(attributes, "(?<=;gene=)[^;]+")]
cds_gff[, gene_id_fbgn := str_extract(attributes, "FBgn[0-9]+")]
cds_gff[, .N, by = .(chr_dm6, strand)]

# x prefix = cds_gff
# i prefix = shared_dt
mel_cds <- cds_gff[
    shared_dt[!is.na(variant.id_mel) & af_mel >= 0.05 & af_mel <= 0.95], 
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

        af_mel = i.af_mel,
        # maf_mel = i.maf_mel,
        effect_mel = i.effect_mel,
        codon_change_mel = i.codon_change_mel,
        aa_change_mel = i.aa_change_mel,
        aa_pos_mel = i.aa_pos_mel,

        variant.id_sim = i.variant.id_sim,
        ref_sim = i.ref_sim,
        alt_sim = i.alt_sim,
        af_sim = i.af_sim,
        # maf_sim = i.maf_sim,
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
        gene_id_fbgn = i.gene_id_fbgn
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
sum(is.na(sim_cds$codon_from_genome)) # should be 0
sum(is.na(sim_cds$codon_ref_annotation)) # should be 0
table(mel_cds$codon_from_genome == mel_cds$codon_ref_annotation)


# codon_genomic_start = genomic coord of 1st base on reference strand
#              _end = 3rd base 
# codon_genomic_raw = the actual codon (eg: CCT) read left to right on reference strand (+)
# codon_from_genome = genomic codon converted into 5' to 3' (reverse complmented for - strand)
# codon_ref_annotation = reference codon annotated by snpeff (from codon_change_sim)

### dm6 codon id: chr:[start]:[end]
mel_cds[, codon_id_dm6_mel := paste(
    chr,
    codon_genomic_start,
    codon_genomic_end,
    sep = ":"
)]

mel_cds[, n_snps_in_codon := .N, by = codon_id_dm6_mel]
mel_1snp_perCodon <- mel_cds[n_snps_in_codon ==1]

shared_dt <- merge(shared_dt, mel_cds[, .(chr, pos, gene_id_fbgn_fromGFF, codon_id_dm6_mel, strand, snp_pos_in_codon)], by=c("chr", "pos"), all.x=T)
setnames(shared_dt, c("strand", "snp_pos_in_codon"), c("strand_mel", "snp_pos_in_codon_mel"))


###############################################################################################
### 3A: D. simulans: Annotate shared table for 1 snp per codon
###############################################################################################
# this annotates all rows TRUE or FALSE depending on if they are in the sim dt 
    # so FALSE includes rows that didn't pass multiple snp per codon filter AND 
    #   those that have no sim annotation at that row
# number of TRUE rows == number of rows in sim_1snp_perCodon

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


# shared_dt[, keep_mel := as.character(variant.id_mel %in% mel_1snp_perCodon$variant.id_mel)]

# shared_dt[, keep_sim := fifelse(
#     is.na(variant.id_sim),
#     NA_character_,
#     as.character(variant.id_sim %in% sim_1snp_perCodon$variant.id_sim)
# )]

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

saveRDS(shared_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.newAF.0.05MAFfilter.annotatedSim.annotatedMel.rds")

###############################################################################################
### 8: Clean up dt
###############################################################################################
filtered_dt <- shared_dt[!is.na(variant.id_mel) | !is.na(variant.id_sim)]
filtered_dt[, maf := NULL]
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
    "n_samps_mel", 
    "af_mel", "maf_mel", 
    "effect_mel", 
    "strand_mel"
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
    "n_samps_sim", 
    "af_sim", "maf_sim", 
    "effect_sim", 
    "flip", "swap", 
    "strand_sim"
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
    #
    "n_samps_mel", "n_samps_sim", 
    "af_mel", "maf_mel", 
    "af_sim", "maf_sim", 
    "effect_mel", "effect_sim",
    "flip", "swap", "strand_sim", "strand_mel", 
    "keep_mel", "keep_sim"

)

new_order <- intersect(new_order, names(filtered_dt))
filtered_dt <- filtered_dt[, ..new_order]
saveRDS(filtered_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.newAF.0.05MAFfilter.annotatedSim.annotatedMel.cleaner.rds")




### inspecting # snps per codons

mel_pos <- filtered_dt[
    !is.na(variant.id_mel),
    .(codon_id = codon_id_dm6_mel,
      mel_pos = snp_pos_in_codon_mel)
]

sim_pos <- filtered_dt[
    !is.na(variant.id_sim),
    .(codon_id = codon_id_dm6_sim,
      sim_pos = snp_pos_in_codon_sim)
]

shared <- merge(mel_pos, sim_pos, by = "codon_id")

shared[, same_position := mel_pos == sim_pos]

table(shared$same_position,useNA="ifany")

shared[, tsp_samePos := same_position]

tsp_map <- unique(shared[, .(
    codon_id = codon_id,
    tsp_samePos
)])

filtered_dt[
    tsp_map,
    on = .(
        codon_id_dm6_mel = codon_id
    ),
    tsp_samePos := i.tsp_samePos
]

filtered_dt[
    tsp_map,
    on = .(
        codon_id_dm6_sim = codon_id
    ),
    tsp_samePos := i.tsp_samePos
]

table(filtered_dt$tsp_samePos, useNA="ifany") 
# TRUE = number of variants with same position variant
# FALSE = 2x number of codons with diff position snps
# NA = species specific variants OR codon id did not match (reading frame was off)


saveRDS(filtered_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.cleaner.1snpCodon.rds")


















### added from pre-classification script. integrate into existing pre-classification script

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

# opp_strand <- filtered_dt[
#     keep_mel=="TRUE" & keep_sim=="TRUE" &
#     !is.na(codon_id_dm6_mel) &
#     !is.na(codon_id_dm6_sim) &
#     !is.na(strand_mel) &
#     !is.na(strand_sim) &
#     strand_mel != strand_sim
# ]
# filtered_dt[codon_id_dm6_mel%in%opp_strand$codon_id_dm6_mel, `:=` (keep_mel="oppStrnd", keep_sim="oppStrnd")]


filtered_dt[keep_mel=="OVERLAP", samePos := "OVERLAP"]
filtered_dt[keep_sim=="OVERLAP", samePos := "OVERLAP"]


### mel cols:
mel_cols <- c(
    "variant.id_mel",
    "ref_mel", "alt_mel", 
    "codon_id_dm6_mel", 
    "codon_ref_mel", "codon_alt_mel",
    "snp_pos_in_codon_mel", 
    "aa_pos_mel", 
    "aa_ref_mel", "aa_alt_mel", 
    "gene_mel", "gene_id_fbgn_fromGFF", 
    "transcript_id_mel", 
    "n_samps_mel", 
    "af_mel", #"maf_mel", 
    "effect_mel", 
    "strand_mel"
)
### sim cols:
sim_cols <- c(
    "variant.id_sim",
    "chr_dsim3", "pos_dsim3", 
    "ref_sim_dm6", "alt_sim_dm6", "ref_sim_dsim3", "alt_sim_dsim3",
    "codon_id_dm6_sim", "codon_id_dsim3",
    "codon_ref_sim", "codon_alt_sim",
    "snp_pos_in_codon_sim", 
    "aa_pos_sim",
    "aa_ref_sim", "aa_alt_sim",
    "gene_sim",
    "transcript_id_sim",
    "n_samps_sim", 
    "af_sim", # "maf_sim",
    "effect_sim", 
    "flip", "swap", 
    "strand_sim"
)
filtered_dt[keep_mel != "TRUE", (mel_cols) := NA]
filtered_dt[keep_sim != "TRUE", (sim_cols) := NA]

# filtered_dt[keep_mel != "TRUE", variant.id_mel := NA_character_]
new_dt <- filtered_dt[!is.na(variant.id_mel) | !is.na(variant.id_sim)]
new_dt <- new_dt[keep_mel == "TRUE" | is.na(keep_mel)]
new_dt <- new_dt[keep_sim == "TRUE" | is.na(keep_sim)]


new_dt[, tsp_samePos := snp_pos_in_codon_mel == snp_pos_in_codon_sim]



new_dt[, .N, by = codon_id_dm6][, table(N)]
new_dt[is.na(tsp_samePos), .N, by = codon_id_dm6][, table(N)] # mixture of diffPos snps and species-specific sites 
new_dt[!is.na(tsp_samePos), .N, by = codon_id_dm6][, table(N)] # all the samePos sites (all TRUE), so all 1 codon only... good!

# View(new_dt[is.na(tsp_samePos), .N, by = codon_id_dm6])
# # get the list of these codons where only 1 codon id dm6 but is not same pos and remove... 

# remove_codons <- new_dt[
#     is.na(tsp_samePos),
#     .N,
#     by = codon_id_dm6
# ][N == 1, codon_id_dm6]

# length(remove_codons)
# new_dt <- new_dt[!(codon_id_dm6 %in% remove_codons)]
### then classify!!!

tsp_map <- new_dt[
    ,
    .(
        pos_mel = snp_pos_in_codon_mel[!is.na(snp_pos_in_codon_mel)][1],
        pos_sim = snp_pos_in_codon_sim[!is.na(snp_pos_in_codon_sim)][1]
    ),
    by = codon_id_dm6
]

tsp_map[, samePos := pos_mel == pos_sim]

new_dt[
    tsp_map,
    on = "codon_id_dm6",
    tsp_samePos := i.samePos
]

# check everything:
table(tsp_map$samePos, useNA = "ifany")
table(new_dt$tsp_samePos, useNA="ifany")


new_dt[
     ,
     .(
         N_rows = .N,
         n_mel = sum(!is.na(variant.id_mel)),
         n_sim = sum(!is.na(variant.id_sim))
     ),
     by = .(codon_id_dm6, tsp_samePos)
 ][
     ,
     .N,
     by = .(tsp_samePos, N_rows, n_mel, n_sim)
 ][order(tsp_samePos, N_rows, n_mel, n_sim)]


 






































 ########################################################
 ### TRASH
 
mel_dt <- readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/DGRP2.source_BCM-HGSC.dm6.snp_dt_SynMissense.10_2026.rds")

mel_gff <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff"


gff <- fread(
     cmd = paste("grep -v '^#'", shQuote(mel_gff)), 
    #  sim_gff, 
     sep="\t", header=F, 
    #  skip = "#",, 
    col.names=c(
        "chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"
    ))

# map chr names to dm6 normal names 
chr_map <- gff[ feature == "region" & grepl("genome=chromosome", attributes),
    .(
        chr_gff = chr,
        chr = sub(".*Name=([^;]+);.*", "\\1", attributes)
    )
]

gff[
     chr_map,
     chr_dm6 := i.chr,
     on = .(chr = chr_gff)
 ]
    
unique(gff[!is.na(chr_dm6), .(chr, chr_dm6)]) # X, 2L, 2R, 3L, 3R, 4, Y

cds_gff <- gff[
    feature == "CDS" & !is.na(chr_dm6), 
    .(chr_dm6, start, end, strand, phase, attributes)
]

# extract genes and transcripts
cds_gff[, transcript_id := str_extract(
    attributes,
    "(?<=Parent=rna-)[^;,]+"
)]

mrna_map <- gff[
    feature == "mRNA",
    .(
        transcript_id = str_extract(attributes, "(?<=ID=rna-)[^;]+"),
        transcript_id_fbtr = str_extract(attributes, "FBtr[0-9]+")
    )
]
cds_gff[
    mrna_map,
    on = "transcript_id",
    transcript_id_fbtr := i.transcript_id_fbtr
]

# cds_gff[, transcript_id_fbtr := str_extract(attributes, "FBtr[0-9]+")]
cds_gff[, gene_name := str_extract(attributes, "(?<=;gene=)[^;]+")]
cds_gff[, gene_id_fbgn := str_extract(attributes, "FBgn[0-9]+")]
cds_gff[, .N, by = .(chr_dm6, strand)]


# names(cds_gff)
# [1] "chr_dm6"       "start"         "end"           "strand"        "phase"        
# [6] "attributes"    "transcript_id" "gene_name"     "gene_id_fbgn" 
# > names(shared_dt)
#  [1] "chr"                  "pos"                  "variant.id_mel"      
#  [4] "ref_mel"              "alt_mel"              "af_mel"              
#  [7] "maf_mel"              "n_samps_mel"          "effect_mel"          
# [10] "impact_mel"           "functional_class_mel" "codon_change_mel"    
# [13] "aa_change_mel"        "aa_length_mel"        "gene_mel"            
# [16] "biotype_mel"          "gene_coding_mel"      "transcript_id_mel"   
# [19] "exon_rank_mel"        "genotype_mel"         "variant.id_sim"      
# [22] "ref_sim"              "alt_sim"              "af_sim"              
# [25] "maf_sim"              "n_samps_sim"          "chr_src"             
# [28] "pos_src"              "ref_src"              "alt_src"             
# [31] "flip"                 "swap"                 "effect_sim"          
# [34] "impact_sim"           "functional_class_sim" "codon_change_sim"    
# [37] "aa_change_sim"        "aa_length_sim"        "gene_sim"            
# [40] "biotype_sim"          "gene_coding_sim"      "transcript_id_sim"   
# [43] "exon_rank_sim"        "genotype_sim"         "gene_id_fbgn"        
# [46] "aa_ref_mel"           "aa_alt_mel"           "aa_ref_sim"          
# [49] "aa_alt_sim"           "aa_pos_mel"           "aa_pos_sim"          
# [52] "codon_ref_mel"        "codon_alt_mel"        "codon_ref_sim"       
# [55] "codon_alt_sim"        "dist_prev_mel"        "dist_next_mel"       
# [58] "dist_prev_sim"        "dist_next_sim"        "keep_sim"  

# x = mel_cds
# i = shared_dt
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

        af_mel = i.af_mel,
        # maf_mel = i.maf_mel,
        effect_mel = i.effect_mel,
        codon_change_mel = i.codon_change_mel,
        aa_change_mel = i.aa_change_mel,
        aa_pos_mel = i.aa_pos_mel,

        variant.id_sim = i.variant.id_sim,
        ref_sim = i.ref_sim,
        alt_sim = i.alt_sim,
        af_sim = i.af_sim,
        # maf_sim = i.maf_sim,
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
        gene_id_fbgn = i.gene_id_fbgn
    )
]

mel_cds[, nt_pos_mel := as.integer(sub(".*c\\.([0-9]+).*", "\\1", aa_change_mel))]
mel_cds[, pos_in_codon := ((nt_pos_mel - 1) %% 3) + 1]
# confirm that this agrees with the snpEff uppercase annotation:
mel_cds[, codon_start_nt := ((nt_pos_mel - 1) %/% 3) * 3 + 1]
mel_cds[, codon_end_nt   := codon_start_nt + 2]

mel_cds[, snp_pos_in_codon := regexpr("[A-Z]", sub("/.*", "", codon_change_mel))]
nrow(mel_cds[pos_in_codon != snp_pos_in_codon]) # should be 0 . good!

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

### check against reference fna


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
nrow(mel_cds[codon_from_genome != codon_ref_annotation])

table(mel_cds$codon_from_genome == mel_cds$codon_ref_annotation)

# codon_genomic_start = genomic coord of 1st base on reference strand
#              _end = 3rd base 
# codon_genomic_raw = the actual codon (eg: CCT) read left to right on reference strand (+)
# codon_from_genome = genomic codon converted into 5' to 3' (reverse complmented for - strand)
# codon_ref_annotation = reference codon annotated by snpeff (from codon_change_sim)


### verified that codon calculated from reference matches that calculated from snp table



### next: create a codon id of sim-based coordinates: chr:[start]:[end]
mel_cds[, codon_id := paste(
    chr,
    codon_genomic_start,
    codon_genomic_end,
    sep = ":"
)]

mel_cds[, n_snps_in_codon := .N, by = codon_id]

mel_1snp_perCodon <- mel_cds[n_snps_in_codon ==1]



### repeat for mel: 1 snp per codon

# saveRDS(mel_cds, "/scratch/ejy4bu/drosophila/inbred/snpDT/DGRP2.source_BCM-HGSC.dm6.snp_dt_SynMissense.codonMap.10_2026.rds")
# saveRDS(mel_1snp_perCodon, "/scratch/ejy4bu/drosophila/inbred/snpDT/DGRP2.source_BCM-HGSC.dm6.snp_dt_SynMissense.1snp_perCodon.rds")

# this annotates all rows TRUE or FALSE depending on if they are in the sim dt 
    # so FALSE includes rows that didn't pass multiple snp per codon filter AND 
    #   those that have no sim annotation at that row
# the number of TRUE rows == number of rows in sim_1snp_perCodon



# there are 6197 rows that are missing from mel cds (annotated as FALSE but have variant id mel)

saveRDS(shared_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.rds")



############################################################################################
### check if any of the sim liftover results show up as in the wrong spot for mel 
    # in order to determine same / different codon snps

shared_dt <- readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.rds")
mel_1snp_perCodon <- readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/DGRP2.source_BCM-HGSC.dm6.snp_dt_SynMissense.1snp_perCodon.rds")
sim_1snp_perCodon <- readRDS(paste0(out_dir, "sim_1snp_perCodon.rds"))

### add dm6 codon id for sim variants:
sim_1snp_perCodon[
    strand == "+",
    `:=`(
        dm6_codon_start = pos - (snp_pos_in_codon - 1),
        dm6_codon_end   = pos + (3 - snp_pos_in_codon)
    )
]

sim_1snp_perCodon[
    strand == "-",
    `:=`(
        dm6_codon_start = pos - (3 - snp_pos_in_codon),
        dm6_codon_end   = pos + (snp_pos_in_codon - 1)
    )
]

sim_1snp_perCodon[
    ,
    dm6_codon_id := paste(
        chr,
        dm6_codon_start,
        dm6_codon_end,
        sep = ":"
    )
]

### rename mel codon id for dm6 coordinates
mel_1snp_perCodon[
    ,
    dm6_codon_id := codon_id
]

### ID codons where 2 sim dsim3 codons map to 1 dm6 codon
dup_dm6_codons <- sim_1snp_perCodon[
    ,
    .N,
    by = dm6_codon_id
][N > 1]


shared_dt[
    variant.id_sim %in% sim_1snp_perCodon[
        dm6_codon_id %in% dup_dm6_codons$dm6_codon_id,
        variant.id_sim
    ],
    keep_sim := "MAP2TO1"
]

table(shared_dt$keep_sim, useNA = "ifany")
shared_dt_all <- shared_dt[
    sim_1snp_perCodon,
    on = .(variant.id_sim = variant.id_sim),
    `:=`(
        strand_sim = i.strand,
        snp_pos_in_codon_sim = i.snp_pos_in_codon,
        codon_id_dsim3 = i.codon_id,
        codon_id_dm6_sim = i.dm6_codon_id
    )
]

shared_dt_all <- shared_dt_all[
    mel_1snp_perCodon,
    on = .(variant.id_mel = variant.id_mel),
    `:=`(
        strand_mel = i.strand,
        snp_pos_in_codon_mel = i.snp_pos_in_codon,
        codon_id_dm6_mel = i.dm6_codon_id,
        gene_id_fbgn_fromGFF=i.gene_id_fbgn_fromGFF, 
        gene_id_fbgn=i.gene_id_fbgn
    )
]
    
# (, .(strand, snp_pos_in_codon, codon_id, dm6_codon_id),)]
shared_dt_all[, classification := NA_character_] # add empty classification colun





shared_dt_cleaner <- shared_dt[, .(
    ### meta stuff 
    variant.id_mel, variant.id_sim,
    chr_dm6=chr, pos_dm6=pos,
    chr_dsim3=chr_src, pos_dsim3=pos_src,
    
    ### nt stuff:
    ref_mel, alt_mel, 
    ref_sim_dm6=ref_sim, alt_sim_dm6=alt_sim, 
    ref_sim_dsim3=ref_src, alt_sim_dsim3=alt_src, 

    ### codon stuff:
    codon_id_dm6_mel, codon_id_dm6_sim, codon_id_dsim3,
    codon_ref_mel, codon_alt_mel, codon_ref_sim, codon_alt_sim,
    snp_pos_in_codon_mel, snp_pos_in_codon_sim,

    ### amino acid stuff:
    aa_pos_mel, aa_pos_sim,
    aa_ref_mel, aa_alt_mel, aa_ref_sim, aa_alt_sim,

    ### gene / transcript stuff:
    gene_mel, gene_sim, 
    gene_id_fbgn,
    gene_id_fbgn_fromGFF,

    transcript_id_mel, transcript_id_sim,

    ### frequency stuff:
    n_samps_mel, n_samps_sim, 
    af_mel, af_sim,
    effect_mel, effect_sim, 

    ### extras:
    flip, swap, 
    strand_mel, strand_sim,
    keep_mel, keep_sim
)]


### clean up rows where keep_mel / sim != TRUE

filtered_dt <- shared_dt_cleaner
filtered_dt[keep_mel != "TRUE", variant.id_mel := NA_integer_]
filtered_dt[keep_sim != "TRUE", variant.id_sim := NA_integer_]

filtered_dt <- filtered_dt[!is.na(variant.id_mel) | !is.na(variant.id_sim)]
### mel cols:
mel_cols <- c(
    "ref_mel", "alt_mel", 
    "codon_id_dm6_mel", 
    "codon_ref_mel", "codon_alt_mel",
    "snp_pos_in_codon_mel", 
    "aa_pos_mel", 
    "aa_ref_mel", "aa_alt_mel", 
    "gene_mel", "gene_id_fbgn_fromGFF", 
    "transcript_id_mel", 
    "n_samps_mel", 
    "af_mel",
    "effect_mel", 
    "strand_mel"
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
    "n_samps_sim", 
    "af_sim", 
    "effect_sim", 
    "flip", "swap", 
    "strand_sim"
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

saveRDS(filtered_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.cleaner.1snpCodon.rds")


### working off filtered_dt


### checking if there are codons that overlap (different reading frame)

