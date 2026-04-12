######### builds a subset RDS of all variants with:
### A) shared polymorphism at same position
### B) multiple polymorphisms within the same codon

#### edit - this is making an RDS heavily favored towards mel. skipping for now

library(data.table)


in_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/all_variants_adjacent_annotation_test.rds"
shared_dt <- readRDS(in_rds)

filtered_dt <- shared_dt[
    (!is.na(ref_mel) & !is.na(ref_sim)) |
    (!is.na(adjacent_var_pos1))
]

message(nrow(shared_dt[!is.na(ref_mel) & !is.na(ref_sim)]), " rows with same site polymorphism")
message(nrow(shared_dt[!is.na(adjacent_var_pos1)]), " rows with same codon polymorphism")
message(nrow(filtered_dt), " total variants kept")

out_csv <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variants_perCodon_sameSite_test500.csv"
subset_table <- filtered_dt[1:500, ]
fwrite(subset_table, out_csv)
message("saved first 500 rows to csv at ", out_csv)



out_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variants_perCodon_sameSite.rds"
fwrite(shared_dt, out_rds)