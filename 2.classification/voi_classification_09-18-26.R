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
# now filtered_dt contains only codons where there is 1 (valid) snp in both species 

# 2 codons have opposite strandedness:
#   3L:18870748:18870750
#   3R:21866692:21866694

filtered_dt[codon_id_dm6_mel=="3L:18870748:18870750" | codon_id_dm6_sim=="3L:18870748:18870750", 
    keep_mel := "oppStrnd"
]
filtered_dt[codon_id_dm6_mel=="3L:18870748:18870750" | codon_id_dm6_sim=="3L:18870748:18870750", 
    keep_sim := "oppStrnd"
]
filtered_dt[codon_id_dm6_mel=="3R:21866692:21866694" | codon_id_dm6_sim=="3R:21866692:21866694", 
    keep_mel := "oppStrnd"
]
filtered_dt[codon_id_dm6_mel=="3R:21866692:21866694" | codon_id_dm6_sim=="3R:21866692:21866694", 
    keep_sim := "oppStrnd"
]


filtered_dt[keep_mel=="OVERLAP", samePos := "OVERLAP"]
filtered_dt[keep_sim=="OVERLAP", samePos := "OVERLAP"]

copy_filtered_dt <- filtered_dt # safe copy 

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
    "af_mel", "maf_mel", 
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
    "af_sim", "maf_sim",
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


new_dt[, codon_id_dm6  := fifelse(
    !is.na(codon_id_dm6_mel),
    codon_id_dm6_mel,
    codon_id_dm6_sim)]

new_dt[, .N, by = codon_id_dm6][, table(N)]
new_dt[is.na(tsp_samePos), .N, by = codon_id_dm6][, table(N)]
new_dt[!is.na(tsp_samePos), .N, by = codon_id_dm6][, table(N)]

View(new_dt[is.na(tsp_samePos), .N, by = codon_id_dm6])
# get the list of these codons where only 1 codon id dm6 but is not same pos and remove... 
### then classify!!!

# filtered_dt[, samePos := as.character(tsp_samePos)]

# filtered_dt[, tsp_samePos := NULL] 



# table(shared$same_position,useNA="ifany")






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
class_dt <- class_dt[keep_mel=="TRUE" | keep_sim == "TRUE"]

class_dt[, codon_id_dm6  := fifelse(
    !is.na(codon_id_dm6_mel),
    codon_id_dm6_mel,
    codon_id_dm6_sim)]

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

class_dt[, class_key := paste(
    samePos,
    total_codons,
    shared_codons,
    total_aa,
    shared_aa,
    mel_aa,
    sim_aa,
    sep = "_"
)]
class_map <- c(

    # same site
    "TRUE_2_2_2_2_2_2" = "A",
    "TRUE_2_2_1_1_1_1" = "B",

    "TRUE_3_1_3_1_2_2" = "C",
    "TRUE_3_1_2_1_1_2" = "D",
    "TRUE_3_1_2_1_2_1" = "E",
    "TRUE_3_1_2_2_2_2" = "F",
    "TRUE_3_1_1_1_1_1" = "G",

    "TRUE_4_0_4_0_2_2" = "H",
    "TRUE_4_0_3_0_1_2" = "I",
    "TRUE_4_0_3_0_2_1" = "J",
    "TRUE_4_0_3_1_2_2" = "K",
    "TRUE_4_0_2_0_1_1" = "L",
    "TRUE_4_0_2_1_1_2" = "M",
    "TRUE_4_0_2_1_2_1" = "N",
    "TRUE_4_0_2_2_2_2" = "O",
    "TRUE_4_0_1_1_1_1" = "P",

    # different site
    "FALSE_4_0_4_0_2_2" = "Q",
    "FALSE_4_0_3_0_1_2" = "R",
    "FALSE_4_0_3_0_2_1" = "S",
    "FALSE_4_0_3_1_2_2" = "T",
    "FALSE_4_0_2_0_1_1" = "U",
    "FALSE_4_0_2_1_1_2" = "V",
    "FALSE_4_0_2_1_2_1" = "W",
    "FALSE_4_0_2_2_2_2" = "X",
    "FALSE_4_0_1_1_1_1" = "Y"
)

same_site <- class_dt[samePos == "TRUE"]
diff_site <- class_dt[samePos == "FALSE"]
diff_site[, dm6_codon_id := fifelse(
    !is.na(codon_id_dm6_mel),
    codon_id_dm6_mel,
    codon_id_dm6_sim
)]
diff_site[, .N, by = codon_id_dm6][, table(N)]

### same site variants: samePos=="TRUE"


class_dt[, class := unname(class_map[class_key])]
same_site[, classification := unname(class_map[class_key])]

table(class_dt$class, useNA = "ifany")
class_dt[is.na(class), 
         .N, 
         by = .(samePos, total_codons, shared_codons,
                total_aa, shared_aa, mel_aa, sim_aa)]
                
class_dt[, .N, by = .(
        class,
        samePos,
        total_codons,
        shared_codons,
        total_aa,
        shared_aa,
        mel_aa,
        sim_aa
    )
][order(class)]

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

