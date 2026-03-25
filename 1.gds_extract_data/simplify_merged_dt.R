library(data.table)
library(dplyr)

out_dir <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_datatables/"

in_rds <- paste0(out_dir, "all_variants_unfiltered.rds")
out_rds <- paste0(out_dir, "all_variants_clean.rds")
out_csv <- paste0(out_dir, "all_variants_clean_500test.csv")
if(!file.exists(out_rds)) file.create(out_rds)
if(!file.exists(out_csv)) file.create(out_csv)

unfiltered_dt <- readRDS(in_rds)

filtered_dt <- unfiltered_dt %>% select(
    chr, pos, 
    ref_mel, alt_mel, ref_sim, alt_sim, 
    af_mel, af_sim, 
    effect_mel, effect_sim, 
    gene_mel, gene_sim, 
    gene_id_mel, gene_id_sim, 
    aa_change_mel, aa_change_sim
)

setDT(filtered_dt)

# extract amino acids from string formatted like "p.Gly->Ala"
filtered_dt[, aa_ref_mel := gsub("^p\\.([[:alpha:]]{3}).*", "\\1" , aa_sub_mel)]
filtered_dt[, aa_alt_mel := gsub("^p\\.[[:alpha:]]{3}->([[:alpha:]]{3}).*", "\\1" , aa_sub_mel)]
filtered_dt[, aa_ref_sim := gsub("^p\\.([[:alpha:]]{3}).*", "\\1" , aa_sub_sim)]
filtered_dt[, aa_alt_sim := gsub("^p\\.[[:alpha:]]{3}->([[:alpha:]]{3}).*", "\\1" , aa_sub_sim)]

filtered_dt[, aa_pos_mel := sub(".*?([0-9]+).*", "\\1", aa_change_mel)]
filtered_dt[, aa_pos_sim := sub(".*?([0-9]+).*", "\\1", aa_change_sim)]

filtered_dt[, c("aa_sub_mel", "aa_sub_sim") := NULL] # remove this column

filtered_dt[, classification := NA_character_] # add empty classification colun

# reorder 
setcolorder(
  filtered_dt,
  c(
    "chr", "pos",
    "ref_mel", "alt_mel", "ref_sim", "alt_sim",
    "aa_ref_mel", "aa_alt_mel", "aa_ref_sim", "aa_alt_sim",
    "classification",
    "af_mel", "af_sim",
    "effect_mel", "effect_sim",
    "gene_mel", "gene_sim",
    "gene_id_mel", "gene_id_sim"
  )
)
# check column names
names(filtered_dt) 

saveRDS(filtered_dt, out_rds)
subset_table <- filtered_dt[1:500, ]
fwrite(subset_table, out_csv)
message("saved first 500 rows to csv at ", out_csv)



# ########### extracting number of variants per codon
# filtered_dt[, aa_pos_mel := sub(".*?([0-9]+).*", "\\1", aa_change_mel)]
# filtered_dt[, aa_pos_sim := sub(".*?([0-9]+).*", "\\1", aa_change_sim)]

# mel_counts <- codon_dt[, .(n_variants = .N), by = .(chr, gene_mel, aa_pos_mel)]
# sim_counts <- codon_dt[, .(n_variants = .N), by = .(chr, gene_sim, aa_pos_sim)]

# mel_summary <- mel_counts[, .(mel = .N), by = n_variants]
# sim_summary <- sim_counts[, .(sim = .N), by = n_variants]

# final_table <- merge(
#    mel_summary,
#    sim_summary,
#    by = "n_variants",
#    all = TRUE
# ) column

# # Rename column
# setnames(final_table, "n_variants", "variants_within_codon")
# final_table[is.na(final_table)] <- 0

# final_table <- final_table[
#    CJ(variants_within_codon = 1:3),
#    on = "variants_within_codon"
# ]

# final_table[is.na(final_table)] <- 0
# fwrite(final_table, "/scratch/ejy4bu/drosophila/gds_analysis/snp_datatables/test_files/variants_per_codon.csv")
 