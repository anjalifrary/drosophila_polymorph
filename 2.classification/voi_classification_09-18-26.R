library(SeqArray)
library(data.table)
library(foreach)
library(doMC)

######################################################################

# ### Pool seq files ###

#######################################################################


### INBRED ###

shared_dt <- readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.cleaner.1snpCodon.rds")

### filter for tsp_samePos TRUE or FALSE (ignore NAs)
filtered_dt <- shared_dt[!is.na(tsp_samePos)]

# now filtered_dt contains only codons where there is 1 (valid) snp in both species 

# 2 codons have opposite strandedness:
#   3L:18870748:18870750
#   3R:21866692:21866694

filtered_dt[, samePos := as.character(tsp_samePos)]

filtered_dt[codon_id_dm6_mel=="3L:18870748:18870750" | codon_id_dm6_sim=="3L:18870748:18870750", 
    samePos := "oppStrnd"
]
filtered_dt[codon_id_dm6_mel=="3R:21866692:21866694" | codon_id_dm6_sim=="3R:21866692:21866694", 
    samePos := "oppStrnd"
]


get_pair <- function(ref,alt) {
    if(is.na(ref) | is.na(alt)) return(NA_character_)
    sort(c(ref,alt)) # returns a vector
}

### i verified that the dm6/dsim3 nt alleles were the same for ALL rows
count_shared <- function(mel, sim){
    if (any(is.na(mel)) | any(is.na(sim))){
        return(c(
            total = NA_integer_,
            shared = NA_integer_,
            Nshared = NA_integer_
        ))
    }
    mel <- unique(mel)
    sim <- unique(sim)
    total <- length(union(mel,sim))
    shared <- length(intersect(mel,sim))
    c(
        total = total,
        shared = shared,
        Nshared = total-shared
    )
}

class_dt <- filtered_dt[samePos == "FALSE" | samePos == "TRUE"]

# get codon & amino acid pairs for mel and sim: 

class_dt[, codon_pair_mel := mapply(get_pair,codon_ref_mel, codon_alt_mel, SIMPLIFY=F)]
class_dt[, codon_pair_sim := mapply(get_pair,codon_ref_sim, codon_alt_sim, SIMPLIFY=F)]

class_dt[, aa_pair_mel := mapply(get_pair, aa_ref_mel, aa_alt_mel, SIMPLIFY=F)]
class_dt[, aa_pair_sim := mapply(get_pair, aa_ref_sim, aa_alt_sim, SIMPLIFY=F)]


### get codon & amino acid counts:
class_dt[, c("total_codons", "shared_codons","Nshared_codons") :=
    transpose(
        mapply(count_shared, codon_pair_mel, codon_pair_sim, SIMPLIFY = FALSE)
    )]

# amino acids:
class_dt[, c("total_aa", "shared_aa", "Nshared_aa") :=
    transpose(
        mapply(count_shared, aa_pair_mel, aa_pair_sim, SIMPLIFY = FALSE)
    )]

class_dt[, mel_aa := lengths(lapply(aa_pair_mel, function(x) unique(na.omit(x))))]
class_dt[, sim_aa := lengths(lapply(aa_pair_sim, function(x) unique(na.omit(x))))]

### same site variants: samePos=="TRUE"


### diff site variants: samePos=="FALSE"


# ### keep the old mapping:
# A  = same_pos 1, total_codons 2, shared_codons 2, total_aa 2, shared_aa 2, mel_aa 2, sim_aa 2
# B  = same_pos 1, total_codons 2, shared_codons 2, total_aa 1, shared_aa 1, mel_aa 1, sim_aa 1

# C  = same_pos 1, total_codons 3, shared_codons 1, total_aa 3, shared_aa 1, mel_aa 2, sim_aa 2
# D  = same_pos 1, total_codons 3, shared_codons 1, total_aa 2, shared_aa 1, mel_aa 1, sim_aa 2
# E  = same_pos 1, total_codons 3, shared_codons 1, total_aa 2, shared_aa 1, mel_aa 2, sim_aa 1
# F  = same_pos 1, total_codons 3, shared_codons 1, total_aa 2, shared_aa 2, mel_aa 2, sim_aa 2
# G  = same_pos 1, total_codons 3, shared_codons 1, total_aa 1, shared_aa 1, mel_aa 1, sim_aa 1

# H  = same_pos 1, total_codons 4, shared_codons 0, total_aa 4, shared_aa 0, mel_aa 2, sim_aa 2
# I  = same_pos 1, total_codons 4, shared_codons 0, total_aa 3, shared_aa 0, mel_aa 1, sim_aa 2
# J  = same_pos 1, total_codons 4, shared_codons 0, total_aa 3, shared_aa 0, mel_aa 2, sim_aa 1
# K  = same_pos 1, total_codons 4, shared_codons 0, total_aa 3, shared_aa 1, mel_aa 2, sim_aa 2
# L  = same_pos 1, total_codons 4, shared_codons 0, total_aa 2, shared_aa 0, mel_aa 1, sim_aa 1
# M  = same_pos 1, total_codons 4, shared_codons 0, total_aa 2, shared_aa 1, mel_aa 1, sim_aa 2
# N  = same_pos 1, total_codons 4, shared_codons 0, total_aa 2, shared_aa 1, mel_aa 2, sim_aa 1
# O  = same_pos 1, total_codons 4, shared_codons 0, total_aa 2, shared_aa 2, mel_aa 2, sim_aa 2
# P  = same_pos 1, total_codons 4, shared_codons 0, total_aa 1, shared_aa 1, mel_aa 1, sim_aa 1

# Q  = same_pos 0, total_codons 4, shared_codons 0, total_aa 4, shared_aa 0, mel_aa 2, sim_aa 2
# R  = same_pos 0, total_codons 4, shared_codons 0, total_aa 3, shared_aa 0, mel_aa 1, sim_aa 2
# S  = same_pos 0, total_codons 4, shared_codons 0, total_aa 3, shared_aa 0, mel_aa 2, sim_aa 1
# T  = same_pos 0, total_codons 4, shared_codons 0, total_aa 3, shared_aa 1, mel_aa 2, sim_aa 2
# U  = same_pos 0, total_codons 4, shared_codons 0, total_aa 2, shared_aa 0, mel_aa 1, sim_aa 1
# V  = same_pos 0, total_codons 4, shared_codons 0, total_aa 2, shared_aa 1, mel_aa 1, sim_aa 2
# W  = same_pos 0, total_codons 4, shared_codons 0, total_aa 2, shared_aa 1, mel_aa 2, sim_aa 1
# X  = same_pos 0, total_codons 4, shared_codons 0, total_aa 2, shared_aa 2, mel_aa 2, sim_aa 2
# Y  = same_pos 0, total_codons 4, shared_codons 0, total_aa 1, shared_aa 1, mel_aa 1, sim_aa 1





#### TRASH #### 
out_dir <- "/scratch/ejy4bu/drosophila/inbred/snpDT/"
in_rds <- paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_clean.rds")
background_rds <- paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.background_1polymorphicCodon.rds")
candidate_rds <- paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.candidate_sharedPolyCodonsOnly.rds")

filtered_dt <- readRDS(in_rds)

# testing on mel only , chr 2L
mel_dt_2L <- filtered_dt[
    !is.na(ref_mel) & chr=="2L", .(
        chr, pos, classification, 
        ref_mel, alt_mel, 
        codon_ref_mel, codon_alt_mel, 
        gene_id_fbgn, gene_mel
    )
]

sim_dt_2L <- filtered_dt[
    !is.na(ref_sim) & chr=="2L", .(
        chr, pos, classification, 
        ref_sim, alt_sim, 
        codon_ref_sim, codon_alt_sim,
        gene_sim
    )
]
unfiltered_dt <- readRDS(paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.rds"))

unfiltered_dt[, aa_pos_mel := gsub(".*?([0-9]+).*", "\\1", aa_change_mel)]
unfiltered_dt[, aa_pos_sim := gsub(".*?([0-9]+).*", "\\1", aa_change_sim)]

mel_dt_2L <- merge(mel_dt_2L, unfiltered_dt[, .(chr, pos, transcript_id_mel, aa_pos_mel)], by=c("chr", "pos"))
same_aa_mel <- mel_dt_2L[
    , if (.N>1) .SD,
    by=.(chr, aa_pos_mel, transcript_id_mel)
]
setorder(mel_dt_2L, chr, pos)
mel_dt_2L[, dist_prev := pos - shift(pos), by = chr]
mel_dt_2L[, dist_next := shift(pos, type = "lead") - pos, by = chr]




sim_dt_2L <- merge(sim_dt_2L, unfiltered_dt[, .(chr, pos, transcript_id_sim, aa_pos_sim)], by=c("chr", "pos"))

same_aa_sim <- sim_dt_2L[
    , if (.N>1) .SD,
    by=.(chr, aa_pos_sim, transcript_id_sim)
]

setorder(sim_dt_2L, chr, pos)
sim_dt_2L[, dist_prev := pos - shift(pos), by = chr]
sim_dt_2L[, dist_next := shift(pos, type = "lead") - pos, by = chr]

