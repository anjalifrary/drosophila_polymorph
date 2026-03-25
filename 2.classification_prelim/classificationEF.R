################################################################################################
# for the classification of same site, same nt, same AND different amino acids 
### distinction: checking for the PAIR of nt and PAIR of aa (which is ref/alt does not matter)
################################################################################################

library(data.table)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_datatables/all_variants_clean.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")

### CLASS E: same position, same nt but swapped ref/alt, same amino acid ############################################
shared_dt[
    ref_mel == alt_sim &
    alt_mel == ref_sim & 
    aa_ref_mel == aa_ref_sim & 
    aa_alt_mel == aa_alt_sim,
    classification := "E"
]
shared_classE <- shared_dt[classification == "E"]
message(nrow(shared_classE), " shared variants in class E")

### CLASS F: same position, same nt but swapped ref/alt, diff amino acid ###########################################
shared_dt[
    ref_mel == alt_sim &
    alt_mel == ref_sim & 
    aa_ref_mel == aa_ref_sim & 
    aa_alt_mel != aa_alt_sim,
    classification := "F"
]
shared_classF <- shared_dt[classification == "F"]
message(nrow(shared_classF), " shared variants in class F")
