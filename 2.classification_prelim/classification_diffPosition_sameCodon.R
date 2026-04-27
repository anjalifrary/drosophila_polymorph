################################################################################################
# for the classification of different site but same codon, same / different nucleotides, same / different amino acids 
################################################################################################

library(data.table)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_variantsOfInterest_qualFiltered.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")


shared_dt[
    !is.na(ref_mel) & # start with mel rows, don't need to do sim bc all should reference each other but checking... 
    (
        adjacent_var_pos1 %in% shared_dt[!is.na(ref_sim), pos] |
        adjacent_var_pos2 %in% shared_dt[!is.na(ref_sim), pos] |
        adjacent_var_pos3 %in% shared_dt[!is.na(ref_sim), pos]
    ),
    classification := ifelse(
    is.na(classification),
    "I",
    paste0(classification, ";I")
    )
]
table(shared_dt$classification)
#  A   A;I     B   B;I     C   C;I     D   D;I     E   E;I     F   F;I     G
# 58653 21971   215   138 22483  7529 14424  5705    52    18  2903  1512  2568
#   G;I     H   H;I     I
#  1371   177    73 58620

#      A  A;I;I      B  B;I;I      C  C;I;I      D  D;I;I      E  E;I;I      F
#  58653  21971    215    138  22483   7529  14424   5705     52     18   2903
#  F;I;I      G  G;I;I      H  H;I;I      I
#   1512   2568   1371    177     73 187142
shared_dt[
    !is.na(ref_sim) & # checking sim, should be same
    (
        adjacent_var_pos1 %in% shared_dt[!is.na(ref_mel), pos] |
        adjacent_var_pos2 %in% shared_dt[!is.na(ref_mel), pos] |
        adjacent_var_pos3 %in% shared_dt[!is.na(ref_mel), pos]
    ),
    classification := ifelse(
    is.na(classification),
    "I",
    paste0(classification, ";I")
    )
]

table(shared_dt$classification) # should be same as before

# clean up
shared_dt[
    ,
    classification := gsub("(^I$)|(;I$)", "", classification)
]
table(shared_dt$classification) # should be same as before

unique_codons_I <- unique(
    shared_dt[grepl("(^|;)I($|;)", classification),
              .(gene_id_mel, codon_start_pos)]
)
message(length(unique_codons_I))
