library(data.table)


in_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/all_variants_adjacent_annotation_test.rds"
shared_dt <- readRDS(in_rds)

filtered_dt <- shared_dt[
    (!is.na(ref_mel) & !is.na(ref_sim)) |
    (!is.na(adjacent_var_pos1))
]
