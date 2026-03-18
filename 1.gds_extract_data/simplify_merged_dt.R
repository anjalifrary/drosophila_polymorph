library(data.table)
library(dplyr)

out_dir <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_datatables/"

in_rds <- paste0(out_dir, "all_variants_unfiltered.rds")
out_rds <- paste0(out_dir, "all_variants_filtered.rds")
out_csv <- paste0(out_dir, "all_variants_filtered_500test.csv")
if(!file.exists(out_rds)) file.create(out_rds)

unfiltered_dt <- readRDS(in_rds)

filtered_dt <- unfiltered_dt %>% select(
    chr, pos, 
    ref_mel, alt_mel, ref_sim, alt_sim, 
    af_mel, af_sim, 
    effect_mel, effect_sim, 
    gene_mel, gene_sim, 
    gene_id_mel, gene_id_sim, 
    aa_sub_mel, aa_sub_sim
)

setDT(filtered_dt)

# extract amino acids from string formatted like "p.Gly->Ala"
filtered_dt[, aa_ref_mel := gsub("^p\\.([[:alpha:]]{3}).*", "\\1" , aa_sub_mel)]
filtered_dt[, aa_alt_mel := gsub("^p\\.[[:alpha:]]{3}->([[:alpha:]]{3}).*", "\\1" , aa_sub_mel)]
filtered_dt[, aa_ref_sim := gsub("^p\\.([[:alpha:]]{3}).*", "\\1" , aa_sub_sim)]
filtered_dt[, aa_alt_sim := gsub("^p\\.[[:alpha:]]{3}->([[:alpha:]]{3}).*", "\\1" , aa_sub_sim)]

filtered_dt[, c("aa_sub_mel", "aa_sub_sim") := NULL] # remove this column

filtered_dt[, classification := NA_character_] # add empty classification colun

# reorder 
setcolorder(
  filtered_dt,
  c(
    "chr", "pos",
    "ref_mel", "alt_mel", "ref_sim", "alt_sim",
    "aa_ref_mel", "aa_alt_mel", "aa_ref_sim", "aa_alt_sim",
    "af_mel", "af_sim",
    "effect_mel", "effect_sim",
    "gene_mel", "gene_sim",
    "gene_id_mel", "gene_id_sim",
    "classification"
  )
)

saveRDS(filtered_dt, out_rds)
subset_table <- filtered_dt[1:500, ]
fwrite(subset_table, out_csv)
message("saved first 500 rows to csv at ", out_csv)