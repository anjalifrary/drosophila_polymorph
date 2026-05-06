library(data.table)

in_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_working.rds"
shared_dt <- readRDS(in_rds)

### A ----- same site variants 
same_site <- shared_dt[!is.na(ref_mel) & !is.na(ref_sim)]

### B ----- different site, same codon, polymorphic in BOTH species:
# join on chr + codon_start_pos,,, one row must be mel and one row must be sim

mel_codons <- unique(shared_dt[!is.na(ref_mel) & is.na(ref_sim), .(chr, codon_start_pos)])
sim_codons <- unique(shared_dt[!is.na(ref_sim) & is.na(ref_mel), .(chr, codon_start_pos)])

shared_codons <- merge(mel_codons, sim_codons, by = c("chr", "codon_start_pos"))

# Keep mel and sim rows whose codon appears in shared_codons
# mel_with_sim_in_codon <- mel_only[shared_codons, on = .(chr, codon_start_pos), nomatch = 0L]
# sim_with_mel_in_codon <- sim_only[shared_codons, on = .(chr, codon_start_pos), nomatch = 0L]

filtered_dt <- rbindlist(list(
  same_site,
  shared_dt[shared_codons, on= .(chr, codon_start_pos)]
), use.names = TRUE)

filtered_dt <- unique(filtered_dt)


message(nrow(same_site), " same-site variants (Condition A)")
message(nrow(shared_codons), " codons polymorphic in both species (Condition B)")

# message(nrow(mel_with_sim_in_codon), " mel variants with sim in same codon (Condition B)")
# message(nrow(sim_with_mel_in_codon), " sim variants with mel in same codon (Condition B)")
message(nrow(filtered_dt), " total unique rows kept")


############# recompute adjacent_var_positions for easy lookup #######################

filtered_dt[, `:=`(
  codon_ref_use = fifelse(!is.na(codon_ref_mel), codon_ref_mel, codon_ref_sim),
  nt_ref_use = fifelse(!is.na(ref_mel), ref_mel, ref_sim)
)]

filtered_dt[, strand := fifelse(
  nt_ref_use == toupper(substr(codon_ref_use, regexpr("[A-Z]", codon_ref_use), 1)), "forward", "reverse"
)]

filtered_dt[strand == "forward", codon_start_pos := pos - (regexpr("[A-Z]", codon_ref_use) - 1)]
filtered_dt[strand == "reverse", codon_start_pos := pos + (regexpr("[A-Z]", codon_ref_use) - 1)]

setkey(filtered_dt, chr, pos)


# Recompute adjacent variant positions
filtered_dt[, aa_pos_toUse := ifelse(!is.na(aa_pos_mel), aa_pos_mel, aa_pos_sim)]

codons_mult_variants <- filtered_dt[
  !is.na(aa_pos_toUse),
  .(n_variants = .N, positions = list(pos)),
  by = .(chr, codon_start_pos)
][n_variants > 1]

message("codons with >1 variant: ", nrow(codons_mult_variants))

# Reset adjacent columns
filtered_dt[, adjacent_var_pos1 := NA_integer_]
filtered_dt[, adjacent_var_pos2 := NA_integer_]
filtered_dt[, adjacent_var_pos3 := NA_integer_]

codons_mult_variants[, `:=`(
  adjacent_variant_pos1 = sapply(positions, \(x) x[1]),
  adjacent_variant_pos2 = sapply(positions, \(x) x[2]),
  adjacent_variant_pos3 = sapply(positions, \(x) x[3])
)]

filtered_dt[codons_mult_variants, on = .(chr, codon_start_pos),
  `:=`(
    adjacent_var_pos1 = i.adjacent_variant_pos1,
    adjacent_var_pos2 = i.adjacent_variant_pos2,
    adjacent_var_pos3 = i.adjacent_variant_pos3
  )
]

# Clean up temp columns
filtered_dt[, c("aa_pos_toUse", "codon_ref_use", "nt_ref_use", "strand") := NULL]
codons_mult_variants[, positions := NULL]

message("Adjacent variant positions recomputed on filtered table")

############### save tables: ###########################

out_csv <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_working.csv"
subset_table <- filtered_dt[1:500, ]
fwrite(subset_table, out_csv)
message("saved first 500 rows to csv at ", out_csv)


out_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_working.rds"
saveRDS(filtered_dt, out_rds)

# # copy to project
# cp "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_working.rds" "/project/berglandlab/anjali/drosophila_polymorphism/mel_sim_sharedTables/currentFiles"
# cp "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_working.csv" "/project/berglandlab/anjali/drosophila_polymorphism/mel_sim_sharedTables/currentFiles"

# in_rds <- "/project/berglandlab/anjali/drosophila_polymorphism/mel_sim_sharedTables/subset_variantsOfInterest_qualFiltered.rds"
# table <- readRDS(in_rds)
# subset <- table[1:500, ]
# fwrite(subset, "/project/berglandlab/anjali/drosophila_polymorphism/mel_sim_sharedTables/subset_variantsOfInterest_qualFiltered.csv")