library(data.table)

in_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variants_perCodon_sameSite.rds"
shared_dt <- readRDS(in_rds)

# subset of variants that are sim
sim_dt <- shared_dt[!is.na(ref_sim)]

# subset of variants that are not in sim (usnig ref_sim so that we dont have dupes)
mel_dt <- shared_dt[is.na(ref_sim)]


# gets the pairs of chr and pos for each adjacent variant (in simulans variants only)
adj_dt <- unique(rbindlist(list(
  sim_dt[!is.na(adjacent_var_pos1), .(chr, pos = adjacent_var_pos1)],
  sim_dt[!is.na(adjacent_var_pos2), .(chr, pos = adjacent_var_pos2)],
  sim_dt[!is.na(adjacent_var_pos3), .(chr, pos = adjacent_var_pos3)]
)))

# set keys for joining tables
setkey(mel_dt, chr, pos)
setkey(adj_dt, chr, pos)

# gets rows where the rows match between shared and adj datatables
# adjacent_matches <- shared_dt[adj_dt, on = .(chr, pos), nomatch = 0L]
adjacent_matches <- mel_dt[adj_dt, on = .(chr, pos), nomatch = 0L]

# join the tables
filtered_dt <- unique(rbindlist(list(sim_dt, adjacent_matches), use.names = TRUE))
setkey(filtered_dt, chr, pos)

message(nrow(sim_dt), " simulans variants")
message(nrow(adjacent_matches), " adjacent mel matches")
message(nrow(filtered_dt), " total rows kept")


out_csv <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variantsOfInterest_test500.csv"
subset_table <- filtered_dt[1:500, ]
fwrite(subset_table, out_csv)
message("saved first 500 rows to csv at ", out_csv)


out_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variantsOfInterest.rds"
saveRDS(filtered_dt, out_rds)


# ######## this creates a table of all the simulans variants out of the filtered table of variants of interest
# ### includes simulans variants with a mel at the same position
# ### includes simulans variants with another variant in the same codon
# ### does not include simulans variants that stand alone / discrete variants

# in_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variants_perCodon_sameSite.rds"
# shared_dt <- readRDS(in_rds)

# filtered_dt <- shared_dt[!is.na(ref_sim)] # filter for all the variants with a variant in simulans

# message(nrow(filtered_dt), " simulans variants kept")

# out_csv <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variants_simulans_test500.csv"
# subset_table <- filtered_dt[1:500, ]
# fwrite(subset_table, out_csv)
# message("saved first 500 rows to csv at ", out_csv)


# out_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variants_simulans.rds"
# saveRDS(filtered_dt, out_rds)