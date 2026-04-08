library(data.table)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/all_variants_clean.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")

# get aa_pos for mel or sim depending on which is available: (update: i think aa_pos matches for mel and sim... not confirmed tho)
shared_dt[, aa_pos_toUse := ifelse(!is.na(aa_pos_mel), aa_pos_mel, aa_pos_sim)]
# and gene id
shared_dt[, gene_id_toUse := ifelse(!is.na(gene_id_mel), gene_id_mel, gene_id_sim)]

# note - gene id will not match for mel and sim... this might be a problem
codons_mult_variants <- shared_dt[
    aa_pos_toUse != "" & !is.na(aa_pos_toUse), # valid amino acid positions
    .(n_variants = .N, positions = list(pos)), # counts number of variants in a single codon and records their positions as list
    by = .(chr, gene_id_toUse, aa_pos_toUse) # merged on chr, gene id, and aa pos 
][n_variants > 1] # for all groups with more than 1 variant per codon 

message("codons with >1 variant: ", nrow(codons_mult_variants)) # 1447067 codons 

# assign adjacent variants to respective nucleotide positions (yes, one will map to self)
codons_mult_variants[, ':='(
    adjacent_variant_pos1 = sapply(positions, \(x) x[1]),
    adjacent_variant_pos2 = sapply(positions, \(x) x[2]),
    adjacent_variant_pos3 = sapply(positions, \(x) x[3]) # NA for many variants
)]

codons_mult_variants[, positions := NULL]

# New columns in shared datatable to store nucleotide position of adjacent variants 
shared_dt[, adjacent_var_pos1 := NA_integer_] 
shared_dt[, adjacent_var_pos2 := NA_integer_] 
shared_dt[, adjacent_var_pos3 := NA_integer_] 

shared_dt[codons_mult_variants, on=.(chr, gene_id_toUse, aa_pos_toUse),
          `:=`(adjacent_var_pos1 = i.adjacent_variant_pos1,
               adjacent_var_pos2 = i.adjacent_variant_pos2,
               adjacent_var_pos3 = i.adjacent_variant_pos3)]

shared_dt[, c("aa_pos_toUse", "gene_id_toUse") := NULL]


# get current column order
cols <- names(shared_dt)

# move classification column to column 3
new_order <- append(
  cols[cols != "classification"],   # all columns except classification
  "classification",
  after = 2                         # put it after column 2 → becomes 3rd column
)
setcolorder(shared_dt, new_order)

out_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/all_variants_adjacent_annotation.rds"
saveRDS(shared_dt, out_rds)
message("Saved updated shared_dt with adjacent variant positions")

out_csv <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/all_variants_adjacent_annotation_test500.csv"
subset_table <- shared_dt[1:500, ]
fwrite(subset_table, out_csv)
message("saved first 500 rows to csv at ", out_csv)
