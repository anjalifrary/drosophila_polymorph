library(data.table)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_variantsOfInterest_qualFiltered.rds")
test_rds_out <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_variantsOfInterest_qualFiltered_testing.rds")
shared_dt <- readRDS(test_rds_out)
# shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")

csv_class <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/classification/classificationScheme.csv")
class_dt <- fread(csv_class)
csv_test <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/classification/classificationScheme_testing.csv")

