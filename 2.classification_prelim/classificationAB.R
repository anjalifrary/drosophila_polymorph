library(data.table)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_datatables/all_variants_clean.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")

### CLASS A: same position, same nt, same amino acid ############################################
shared_dt[
    ref_mel == ref_sim &
    alt_mel == alt_sim & 
    aa_ref_mel == aa_ref_sim & 
    aa_alt_mel == aa_alt_sim,
    classification := "A"
]
shared_classA <- shared_dt[classification == "A"]
message(nrow(shared_classA), " shared variants with same pos, nt, and aa")

### CLASS B: same position, same nt, diff amino acid ###########################################
shared_dt[
    ref_mel == ref_sim &
    alt_mel == alt_sim & 
    aa_ref_mel == aa_ref_sim & 
    aa_alt_mel != aa_alt_sim,
    classification := "B"
]

# outputs summary of classification stats 
table(shared_dt$classification, useNA = "ifany")

saveRDS(shared_dt, rds_file)

### updating csv for viewability
subset_table <- shared_dt[1:500, ]
table(subset_table$classification, useNA = "ifany")
fwrite(subset_table, "/scratch/ejy4bu/drosophila/gds_analysis/snp_datatables/all_variants_clean_500test.csv")
