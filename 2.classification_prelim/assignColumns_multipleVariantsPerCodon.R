library(data.table)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/all_variants_clean.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")

# get aa_pos for mel or sim depending on which is available: (update: i think aa_pos matches for mel and sim... not confirmed tho)
shared_dt[, aa_pos_toUse := ifelse(!is.na(aa_pos_mel), aa_pos_mel, aa_pos_sim)]

############# new way based on calculating start nt position for each codon
# message((7684 - 3) %/% 3 * 3 + 3)

# shared_dt[, codon_start_pos := (pos - 3) %/% 3 * 3 + 3]
# shared_dt[, codon_start_pos := pos - (regexpr("[A-Z]", codon_ref_mel) - 1)]

# assign a codon reference to use and a nt reference to use 
shared_dt[, `:=`(
  codon_ref_use = fifelse(!is.na(codon_ref_mel), codon_ref_mel, codon_ref_sim),
  nt_ref_use = fifelse(!is.na(ref_mel), ref_mel, ref_sim)
)]

# determine strand
shared_dt[, strand := fifelse(
  nt_ref_use == toupper(substr(codon_ref_use, regexpr("[A-Z]", codon_ref_use), 1)), "forward", "reverse"
)]

# calculate codon_start_pos for forward and reverse strands
shared_dt[strand == "forward", codon_start_pos := pos - (regexpr("[A-Z]", codon_ref_use) - 1)]
shared_dt[strand == "reverse", codon_start_pos := pos + (regexpr("[A-Z]", codon_ref_use) - 1)]

codons_mult_variants <- shared_dt[
    !is.na(aa_pos_toUse), # valid amino acid positions
    .(n_variants = .N, positions = list(pos)), # counts number of variants in a single codon and records their positions as list
    by = .(chr, codon_start_pos) # merged on chr, and codon start pos
][n_variants > 1] # for all groups with more than 1 variant per codon 

message("codons with >1 variant: ", nrow(codons_mult_variants)) # 779373 codons, now 713565

# New columns in shared datatable to store nucleotide position of adjacent variants 
shared_dt[, adjacent_var_pos1 := NA_integer_] 
shared_dt[, adjacent_var_pos2 := NA_integer_] 
shared_dt[, adjacent_var_pos3 := NA_integer_] 

# assign adjacent variants to respective nucleotide positions (yes, one will map to self)
codons_mult_variants[, ':='(
    adjacent_variant_pos1 = sapply(positions, \(x) x[1]),
    adjacent_variant_pos2 = sapply(positions, \(x) x[2]),
    adjacent_variant_pos3 = sapply(positions, \(x) x[3]) # NA for many variants
)]

shared_dt[codons_mult_variants, on=.(chr, codon_start_pos),
        `:=` (adjacent_var_pos1 = i.adjacent_variant_pos1,
               adjacent_var_pos2 = i.adjacent_variant_pos2,
               adjacent_var_pos3 = i.adjacent_variant_pos3)
]


codons_mult_variants[, positions := NULL]


shared_dt[, c("aa_pos_toUse") := NULL]



################ fix column orders ######################################

# get current column order
cols <- names(shared_dt)

# move classification column to column 3
new_order <- append(
  cols[cols != "classification"],   # all columns except classification
  "classification",
  after = 2                         # put it after column 2 → becomes 3rd column
)
setcolorder(shared_dt, new_order)

cols <- names(shared_dt)
new_order <- append(
  cols[cols != "codon_start_pos"],   # all columns except codon_start
  "codon_start_pos",
  after = 2                         # put it after column 2 → becomes 3rd column
)
setcolorder(shared_dt, new_order)



################## save files ############################################

out_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/all_variants_adjacent_annotation_test.rds"
saveRDS(shared_dt, out_rds)
message("Saved updated shared_dt with adjacent variant positions")

out_csv <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/all_variants_adjacent_annotation_test_test500.csv"
subset_table <- shared_dt[1:500, ]
fwrite(subset_table, out_csv)
message("saved first 500 rows to csv at ", out_csv)


