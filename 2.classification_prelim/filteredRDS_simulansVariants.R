library(data.table)


in_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variants_perCodon_sameSite.rds"
shared_dt <- readRDS(in_rds)

filtered_dt <- shared_dt[!is.na(ref_sim)] # filter for all the variants with a variant in simulans

message(nrow(filtered_dt), " simulans variants kept")

out_csv <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variants_simulans_test500.csv"
subset_table <- filtered_dt[1:500, ]
fwrite(subset_table, out_csv)
message("saved first 500 rows to csv at ", out_csv)


out_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variants_simulans.rds"
saveRDS(filtered_dt, out_rds)