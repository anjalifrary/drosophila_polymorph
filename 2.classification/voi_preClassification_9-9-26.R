library(SeqArray)
library(data.table)
library(foreach)
library(doMC)

######################################################################

# ### Pool seq files ###

#######################################################################

### INBRED ###
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

