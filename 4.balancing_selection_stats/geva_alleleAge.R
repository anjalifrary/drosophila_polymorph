library(data.table)

geva <- fread("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/AlleleAges.VA.cm_GEVA.txt")


rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total rows in shared table")


# extract chr from id 
geva[, chr := tstrsplit(id, ".", fixed = TRUE)[[1]]]

# table(geva$chr, useNA = "ifany")

setnames(geva, "position", "pos")

# merge geva and rds file 

shared_dt <- merge(
    shared_dt, geva[, .(chr, pos, PostMedian, PostMode)],
    by = c("chr", "pos"),
    all.x=T
)

# number of rows that are missing age data.table
sum(!is.na(shared_dt$PostMedian))


out_rds <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed_geva.rds")

saveRDS(shared_dt, out_rds)
message("classified RDS written to: ", out_rds)

# CSV SUBSET
subset_table <- shared_dt[1:500, ]
# table(subset_table$classification, useNA = "ifany")
fwrite(subset_table, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed_age_test500.csv")
