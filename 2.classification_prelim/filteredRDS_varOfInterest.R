library(data.table)

in_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/merged_tables/quality/all_quality_variants_MAF5_clean.rds"
# in_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed_geva_MAF5.rds"
shared_dt <- readRDS(in_rds)
# shared_dt <- voi
### 1. confirm correct codon start pos 

shared_dt[, `:=`(
  codon_ref_use = fifelse(!is.na(codon_ref_mel), codon_ref_mel, codon_ref_sim),
  nt_ref_use = fifelse(!is.na(ref_mel), ref_mel, ref_sim)
)]

shared_dt[, strand := fifelse(
  nt_ref_use == toupper(substr(codon_ref_use, regexpr("[A-Z]", codon_ref_use), 1)), "forward", "reverse"
)]

shared_dt[strand == "forward", codon_start_pos := pos - (regexpr("[A-Z]", codon_ref_use) - 1)]
shared_dt[strand == "reverse", codon_start_pos := pos + (regexpr("[A-Z]", codon_ref_use) - 1)]
shared_dt[, c("codon_ref_use", "nt_ref_use", "strand") := NULL]
# setkey(filtered_dt, chr, pos)

### 2: same-site variants, defined by mel and sim are both present at nucleotide positions
same_site_all <- shared_dt[!is.na(ref_mel) & !is.na(ref_sim)]
# ensure that there is only 1 same-site polymorphism per codon
same_site_codon_counts <- same_site_all[, .(n_same = .N), by = .(chr, codon_start_pos)]

same_site <- same_site_all[
  same_site_codon_counts[n_same == 1],
  on = .(chr, codon_start_pos),
  nomatch = 0L
]

message(nrow(same_site), " same-site variants")

### 3: different-site variants, defined by 1 polymorphic site per species per codon
# get mel and sim unique rows (avoid duplicating same-site variants)
mel_only <- shared_dt[!is.na(ref_mel) & is.na(ref_sim)]
sim_only <- shared_dt[!is.na(ref_sim) & is.na(ref_mel)]

# count: varants per species per codon
mel_codon_count <- mel_only[, .(n_mel = .N), by = .(chr, codon_start_pos)]
sim_codon_count <- sim_only[, .(n_sim = .N), by = .(chr, codon_start_pos)]

# filter for codons where each species is polymorphic at a single site per codon
shared_codons_1 <- merge(mel_codon_count[n_mel == 1], sim_codon_count[n_sim == 1],
  by = c("chr", "codon_start_pos")
)
message(nrow(shared_codons_1), " codons with 1 mel and 1 sim variant")

# merge with the mel only and sim only rows
diff_site <- rbindlist(list(
  mel_only[shared_codons_1, on = .(chr, codon_start_pos), nomatch = 0L],
  sim_only[shared_codons_1, on = .(chr, codon_start_pos), nomatch = 0L]
), use.names = TRUE)
message(nrow(diff_site), " different site variants (2 x # codons)")

diff_site[, c("n_mel", "n_sim") := NULL]

# remove rows where there is a polymorphism at same site and different site within the codon (= all three codon positions)
both_conditions <- fintersect(
  same_site[, .(chr, codon_start_pos)],
  shared_codons_1[, .(chr, codon_start_pos)]
)
message(nrow(both_conditions), " codons appearing in both A and B — dropping")

same_site  <- same_site[!both_conditions,  on = .(chr, codon_start_pos)]
diff_site  <- diff_site[!both_conditions,  on = .(chr, codon_start_pos)]
message(nrow(same_site), " same site variants and ", nrow(diff_site)/2, " diff site codons, (variants)", nrow(diff_site))

same_site[, c("n_same") := NULL]
### 4. combine same and diff site variants
filtered_dt <- rbindlist(list(same_site, diff_site), use.names=TRUE)
message(nrow(filtered_dt), " same & diff site merged rows before calling unique")
filtered_dt <- unique(filtered_dt)
message(nrow(filtered_dt), " after unique, same & diff site rows")

setkey(filtered_dt, chr, pos)


### 5. Recompute adjacent variant positions (remove null pointers)
filtered_dt[, aa_pos_toUse := ifelse(!is.na(aa_pos_mel), aa_pos_mel, aa_pos_sim)]

codons_mult_variants <- filtered_dt[
  !is.na(aa_pos_toUse),
  .(n_variants = .N, positions = list(pos)), #list(sort(pos)) ? 
  by = .(chr, codon_start_pos)
][n_variants > 1]
message("codons with >1 variant in filtered table (should be == diff-site): ", nrow(codons_mult_variants))

filtered_dt[, aa_pos_toUse := NULL]


filtered_dt[, c("adjacent_var_pos1", "adjacent_var_pos2", "adjacent_var_pos3") := NA_integer_]

codons_mult_variants[, `:=`(
  adjacent_variant_pos1 = sapply(positions, \(x) x[1]),
  adjacent_variant_pos2 = sapply(positions, \(x) x[2]),
  adjacent_variant_pos3 = sapply(positions, \(x) x[3])
)]

filtered_dt[codons_mult_variants, on = .(chr, codon_start_pos), `:=`(
  adjacent_var_pos1 = i.adjacent_variant_pos1,
  adjacent_var_pos2 = i.adjacent_variant_pos2,
  adjacent_var_pos3 = i.adjacent_variant_pos3
)]

filtered_dt[, aa_pos_toUse := NULL]
codons_mult_variants[, positions := NULL]

message("Adjacent variant positions recomputed on filtered table")


############### save tables: ###########################

# out_csv <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_final.csv"
# subset_table <- filtered_dt[1:500, ]
# fwrite(subset_table, out_csv)
# message("saved first 500 rows to csv at ", out_csv)


out_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/subset_fromBG_qualVar_ofInterest_MAF5_06-29-2026.rds"
saveRDS(filtered_dt, out_rds)

# cp /scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_final.rds /project/berglandlab/anjali/drosophila_polymorphism/mel_sim_sharedTables