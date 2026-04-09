library(data.table)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/all_variants_clean.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")

