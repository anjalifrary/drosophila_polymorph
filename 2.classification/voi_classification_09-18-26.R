library(SeqArray)
library(data.table)
library(foreach)
library(doMC)

### to run after preClassification


dataset <- "inbred"    # poolseq or inbred

if (dataset == "poolseq") {
    proj_dir <- "/project/berglandlab/anjali/drosophila_polymorphism/classification/poolseq_dest/classed/"
    # input file
    shared_dt <- readRDS("/scratch/ejy4bu/drosophila/DEST_remake/snpDT/shared/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.annotatedSim.annotatedMel.cleaner.rds")

    # output files
    final_dt_save <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.shared.bothMelSim.classed.rds"
    sharedOnly_dt <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.shared.classed.rds"
    candidates_dt <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.candidatesABFGOPXY.classed.rds"
    csv_class <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.classification.csv"

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
        "effect_mel", 
        "strand_mel", 
        "global_af_mel", "poly_af_mel", "poly_maf_mel",
        "nLocales_poly_mel", "nLocales_total_mel", "nLocales_missing_mel", 
        "nLocales_poly_proj_round_mel", "nLocales_poly_proj_binomial_mel"

    )
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
        "effect_sim", 
        "flip", "swap", 
        "strand_sim", 
        "global_af_sim", "poly_af_sim", "poly_maf_sim",
        "nLocales_poly_sim", "nLocales_total_sim", "nLocales_missing_sim", 
        "nLocales_poly_proj_round_sim", "nLocales_poly_proj_binomial_sim"
    )

} else if (dataset == "inbred") {
    proj_dir <- "/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/"

    shared_dt <- readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.newAF.0.05MAFfilter.annotatedSim.annotatedMel.cleaner.rds")
    # shared_dt <- readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.newAF.0.01MAFfilter.annotatedSim.annotatedMel.cleaner.rds")

    final_dt_save <- "/scratch/ejy4bu/drosophila/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.11_2026.shared.bothMelSim.0.05MAFfilter.classed.rds"
    sharedOnly_dt <- "/scratch/ejy4bu/drosophila/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.11_2026.shared.0.05MAFfilter.classed.rds"
    candidates_dt <- "/scratch/ejy4bu/drosophila/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.11_2026.candidatesABFGOPXY.0.05MAFfilter.classed.rds"    
    csv_class <- "/scratch/ejy4bu/drosophila/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.classification.11_2026.newAF.0.05MAF.csv"

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
} else {
    stop("invalid species")
}


filtered_dt <- shared_dt

### look for opposite strand codons:
filtered_dt[, codon_id_dm6  := fifelse(
    !is.na(codon_id_dm6_mel),
    codon_id_dm6_mel,
    codon_id_dm6_sim)]

# realized I filtered for opp codons earlier but missed some... 
opp_codons <- filtered_dt[
    ,
    .(
        has_plus  = any(strand_mel == "+" | strand_sim == "+", na.rm = TRUE),
        has_minus = any(strand_mel == "-" | strand_sim == "-", na.rm = TRUE)
    ),
    by = codon_id_dm6
][has_plus & has_minus, codon_id_dm6]

filtered_dt[codon_id_dm6%in%opp_codons, `:=` (keep_mel = "oppStrand", keep_sim ="oppStrand")]

filtered_dt[keep_mel != "TRUE", (mel_cols) := NA]
filtered_dt[keep_sim != "TRUE", (sim_cols) := NA]

nrow(filtered_dt[!is.na(variant.id_mel) & !is.na(variant.id_sim)])
nrow(filtered_dt[!is.na(variant.id_mel) | !is.na(variant.id_sim)])

# filtered_dt[keep_mel != "TRUE", variant.id_mel := NA_character_]
new_dt <- filtered_dt[!is.na(variant.id_mel) | !is.na(variant.id_sim)]
new_dt <- new_dt[keep_mel == "TRUE" | keep_mel == "No_mel_variant"]
new_dt <- new_dt[keep_sim == "TRUE" | keep_sim == "No_sim_variant"]
nrow(new_dt)

new_dt[, tsp_samePos := snp_pos_in_codon_mel == snp_pos_in_codon_sim]



new_dt[, .N, by = codon_id_dm6][, table(N)]
new_dt[is.na(tsp_samePos), .N, by = codon_id_dm6][, table(N)] # mixture of diffPos snps and species-specific sites 
new_dt[!is.na(tsp_samePos), .N, by = codon_id_dm6][, table(N)] # all the samePos sites (all TRUE), so all 1 codon only... good!

tsp_map <- new_dt[, .(
        pos_mel = snp_pos_in_codon_mel[!is.na(snp_pos_in_codon_mel)][1],
        pos_sim = snp_pos_in_codon_sim[!is.na(snp_pos_in_codon_sim)][1]
    ), by = codon_id_dm6
]

tsp_map[, samePos := pos_mel == pos_sim]

new_dt[ tsp_map, on = "codon_id_dm6", tsp_samePos := i.samePos]

# check everything:
table(tsp_map$samePos, useNA = "ifany")
table(new_dt$tsp_samePos, useNA="ifany")


new_dt[, .(
         N_rows = .N,
         n_mel = sum(!is.na(variant.id_mel)),
         n_sim = sum(!is.na(variant.id_sim))
     ), by = .(codon_id_dm6, tsp_samePos)
 ][, .N, by = .(tsp_samePos, N_rows, n_mel, n_sim)][order(tsp_samePos, N_rows, n_mel, n_sim)]
#  tsp_samePos N_rows n_mel n_sim       N
#         <lgcl>  <int> <int> <int>   <int>
# 1:       FALSE      2     1     1  174353 --- diffPos
# 2:        TRUE      1     1     1  126095 --- samePos
# 3:          NA      1     0     1  276734 --- simOnly
# 4:          NA      1     1     0 2682971 --- melOnly

# # inbred:
# 1:       FALSE      2     1     1  23979
# 2:        TRUE      1     1     1 101938
# 3:          NA      1     0     1 866871
# 4:          NA      1     1     0 344098


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

# new_dt[, samePos := NULL]
setnames(new_dt, "tsp_samePos", "samePos")

class_dt <- new_dt[samePos == "FALSE" | samePos == "TRUE"]
nrow(class_dt)==nrow(new_dt) # should be TRUE or FALSE???? 

class_dt <- class_dt[keep_mel=="TRUE" | keep_sim == "TRUE"]

### retaining the same lettering scheme from the first analysis:
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

# split into same and diff site within codon snps:
same_site <- class_dt[samePos == TRUE]
diff_site <- class_dt[samePos == FALSE]

same_site[, .N, by = codon_id_dm6][, table(N)] # all 1
diff_site[, .N, by = codon_id_dm6][, table(N)] # all 2

### same site classification:
same_site[, codon_pair_mel := mapply(get_pair, codon_ref_mel, codon_alt_mel, SIMPLIFY=F)]
same_site[, codon_pair_sim := mapply(get_pair, codon_ref_sim, codon_alt_sim, SIMPLIFY=F)]

same_site[, aa_pair_mel := mapply(get_pair, aa_ref_mel, aa_alt_mel, SIMPLIFY=F)]
same_site[, aa_pair_sim := mapply(get_pair, aa_ref_sim, aa_alt_sim, SIMPLIFY=F)]

### get codon & amino acid counts:
same_site[, c("total_codons", "shared_codons","Nshared_codons") :=
    transpose(mapply(count_shared, codon_pair_mel, codon_pair_sim, SIMPLIFY = FALSE))]

same_site[, c("total_aa", "shared_aa", "Nshared_aa") :=
    transpose(mapply(count_shared, aa_pair_mel, aa_pair_sim, SIMPLIFY = FALSE))]

same_site[, mel_aa := lengths(lapply(aa_pair_mel, function(x) unique(na.omit(x))))]
same_site[, sim_aa := lengths(lapply(aa_pair_sim, function(x) unique(na.omit(x))))]

same_site[, class_key := paste(
    samePos,
    total_codons,
    shared_codons,
    total_aa,
    shared_aa,
    mel_aa,
    sim_aa,
    sep = "_"
)]
same_site[, classification := unname(class_map[class_key])]

### diff site:
diff_codons <- diff_site[
    ,
    .(
        codon_ref_mel = codon_ref_mel[!is.na(codon_ref_mel)][1],
        codon_alt_mel = codon_alt_mel[!is.na(codon_alt_mel)][1],
        aa_ref_mel    = aa_ref_mel[!is.na(aa_ref_mel)][1],
        aa_alt_mel    = aa_alt_mel[!is.na(aa_alt_mel)][1],

        codon_ref_sim = codon_ref_sim[!is.na(codon_ref_sim)][1],
        codon_alt_sim = codon_alt_sim[!is.na(codon_alt_sim)][1],
        aa_ref_sim    = aa_ref_sim[!is.na(aa_ref_sim)][1],
        aa_alt_sim    = aa_alt_sim[!is.na(aa_alt_sim)][1]
    ),
    by = codon_id_dm6
]
diff_codons[, codon_pair_mel := mapply(get_pair, codon_ref_mel, codon_alt_mel, SIMPLIFY=F)]
diff_codons[, codon_pair_sim := mapply(get_pair, codon_ref_sim, codon_alt_sim, SIMPLIFY=F)]

diff_codons[, aa_pair_mel := mapply(get_pair, aa_ref_mel, aa_alt_mel, SIMPLIFY=F)]
diff_codons[, aa_pair_sim := mapply(get_pair, aa_ref_sim, aa_alt_sim, SIMPLIFY=F)]

### get codon & amino acid counts:
diff_codons[, c("total_codons", "shared_codons","Nshared_codons") :=
    transpose(mapply(count_shared, codon_pair_mel, codon_pair_sim, SIMPLIFY = FALSE))]

diff_codons[, c("total_aa", "shared_aa", "Nshared_aa") :=
    transpose(mapply(count_shared, aa_pair_mel, aa_pair_sim, SIMPLIFY = FALSE))]

diff_codons[, mel_aa := lengths(lapply(aa_pair_mel, function(x) unique(na.omit(x))))]
diff_codons[, sim_aa := lengths(lapply(aa_pair_sim, function(x) unique(na.omit(x))))]

diff_codons[
    ,
    .N,
    by = .(
        total_codons,
        shared_codons,
        Nshared_codons,
        total_aa,
        shared_aa,
        Nshared_aa,
        mel_aa,
        sim_aa
    )
][order(
    total_codons,
    shared_codons,
    total_aa,
    shared_aa,
    mel_aa,
    sim_aa
)]

diff_codons[, class_key := paste(
    FALSE,
    total_codons,
    shared_codons,
    total_aa,
    shared_aa,
    mel_aa,
    sim_aa,
    sep = "_"
)]
diff_codons[, classification := unname(class_map[class_key])]

View(same_site)
View(diff_codons)

### merge same_site and diff_codons back to final_dt
class_lookup <- rbindlist(
    list(
        same_site[, .(codon_id_dm6, classification)],
        diff_codons[, .(codon_id_dm6, classification)]
    ),
    use.names = TRUE
)

# make sure each codon has only one classification
class_lookup[, .N, by = codon_id_dm6][N > 1]

final_dt <- new_dt[
    class_lookup,
    on = "codon_id_dm6",
    classification := i.classification
]
table(final_dt$classification, useNA = "ifany")

final_dt[
    !is.na(classification),
    .N,
    by = classification
][order(classification)]

saveRDS(final_dt, final_dt_save)

sharedOnly <- final_dt[!is.na(classification)]
saveRDS(sharedOnly, sharedOnly_dt)

candidates <- final_dt[classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y")]
saveRDS(candidates, candidates_dt)


class_table <- rbindlist(
    list(
        same_site[, .(
            Classification = classification,
            same_pos = TRUE,
            total_codons,
            shared_codons,
            Nshared_codons,
            total_aa,
            shared_aa,
            Nshared_aa,
            mel_aa,
            sim_aa
        )],
        diff_codons[, .(
            Classification = classification,
            same_pos = FALSE,
            total_codons,
            shared_codons,
            Nshared_codons,
            total_aa,
            shared_aa,
            Nshared_aa,
            mel_aa,
            sim_aa
        )]
    ),
    use.names = TRUE
)

class_table <- class_table[
    ,
    .(
        Count = .N,
        same_pos = same_pos[1],
        total_codons = total_codons[1],
        shared_codons = shared_codons[1],
        Nshared_codons = Nshared_codons[1],
        total_aa = total_aa[1],
        shared_aa = shared_aa[1],
        Nshared_aa = Nshared_aa[1],
        mel_aa = mel_aa[1],
        sim_aa = sim_aa[1]
    ),
    by = Classification
]

setorder(class_table, Classification)

class_table[, `Convergent/Divergent` := fifelse(
    Classification %in% c("A", "B", "F", "G", "O", "P", "X", "Y"),
    "Convergent",
    "Divergent"
)]

class_table[, `Shared/Independent` := fifelse(
    Classification %in% c("A", "B"),
    "Either",
    "Independent"
)]
setcolorder(
    class_table,
    c(
        "Classification",
        "Count",
        "same_pos",
        "total_codons",
        "shared_codons",
        "Nshared_codons",
        "total_aa",
        "shared_aa",
        "Nshared_aa",
        "mel_aa",
        "sim_aa",
        "Convergent/Divergent",
        "Shared/Independent"
    )
)


fwrite(class_table, csv_class)

message("classification table written to: ", csv_class)


file.copy(c(final_dt_save, sharedOnly_dt, candidates_dt, csv_class), proj_dir, overwrite=T)