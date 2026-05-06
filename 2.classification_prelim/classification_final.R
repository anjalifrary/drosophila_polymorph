library(data.table)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_working.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")

csv_class <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/classification/classification_table_current.csv")


# 1. filter for codons where mel and sim are both polymorphic at ONE site within the codon
# this includes same site (same row) and different site (different rows, same codon_start_pos)


# shared_dt <- unique(shared_dt, by = c("chr", "pos"))
# message(nrow(shared_dt))

# 2. same site classification: rows where ref_mel and ref_sim are not empty

# 3. different site classification: 

# 4. write classification table 
# same_pos (1 = same, 0 = different) | total_uniqueCodons | shared_codons | Nshared_codons | total_uniqueAA | shared_AA | Nshared_AA

# save classification table to csv_class
# fill in classification column in rds, save to rds file