library(data.table)
library(dplyr)


# ######################################################################
# # ### Pool seq files ###
# out_dir <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/merged_tables/quality/"
# in_rds <- paste0(out_dir, "all_quality_variants_merge_unfilt.rds")
# out_rds <- paste0(out_dir, "all_quality_variants_clean.rds")
# # out_csv <- paste0(out_dir, "all_quality_variants_clean_500test.csv")
# if(!file.exists(out_rds)) file.create(out_rds)
# # if(!file.exists(out_csv)) file.create(out_csv)

#######################################################################
### INBRED ###
out_dir <- "/scratch/ejy4bu/drosophila/inbred/snpDT/"
unfiltered_dt <- readRDS(paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.rds"))
out_rds <- paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_clean.rds")

library(AnnotationDbi)
library(org.Dm.eg.db)
gene_map <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys = unique(unfiltered_dt$gene_mel),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "FLYBASE")
)

gene_map <- as.data.table(gene_map)
setnames(gene_map, c("SYMBOL", "FLYBASE"), c("gene_mel", "gene_id_fbgn"))

gene_map <- unique(gene_map[, .(gene_mel, gene_id_fbgn)])

unfiltered_dt[gene_map, gene_id_fbgn := i.gene_id_fbgn, on = "gene_mel"]

# names(unfiltered_dt)
#  [1] "gene_mel"             "chr"                  "pos"                 
#  [4] "variant.id_mel"       "ref_mel"              "alt_mel"             
#  [7] "af_mel"               "maf_mel"              "n_samps_mel"         
# [10] "effect_mel"           "impact_mel"           "functional_class_mel"
# [13] "codon_change_mel"     "aa_change_mel"        "aa_length_mel"       
# [16] "biotype_mel"          "gene_coding_mel"      "transcript_id_mel"   
# [19] "exon_rank_mel"        "genotype_mel"         
# "variant.id_sim"      
# [22] "ref_sim"              "alt_sim"              "af_sim"              
# [25] "maf_sim"              "n_samps_sim"          "effect_sim"          
# [28] "impact_sim"           "functional_class_sim" "codon_change_sim"    
# [31] "aa_change_sim"        "aa_length_sim"        "gene_sim"            
# [34] "biotype_sim"          "gene_coding_sim"      "transcript_id_sim"   
# [37] "exon_rank_sim"        "genotype_sim"         "gene_id_fbgn"        
# > 
# don't keep: impact, functional_class, biotype, aa_length, biotype, geen_coding, transcript_id, exon_rank
  # = remove 16 columns (keep 23)

filtered_dt <- unfiltered_dt[, .(
    chr, pos, 
    variant.id_mel, variant.id_sim,
    n_samps_mel, n_samps_sim,
    ref_mel, alt_mel, ref_sim, alt_sim, 
    af_mel, maf_mel, af_sim, maf_sim, 
    effect_mel, effect_sim, 
    gene_id_fbgn, gene_mel, gene_sim, 
    # gene_id_mel, gene_id_sim, 
    # nt_change_mel, nt_change_sim,
    codon_change_mel, codon_change_sim,
    aa_change_mel, aa_change_sim
)]

setDT(filtered_dt)

##### extract amino acids from string formatted like "p.Ser795Ile"
filtered_dt[, aa_ref_mel := gsub("^p\\.([[:alpha:]]{3}).*", "\\1" , aa_change_mel)]   # Ser
# filtered_dt[, aa_alt_mel := gsub(".*([[:alpha:]]{3})$", "\\1" , aa_change_mel)]       # Ile
filtered_dt[, aa_alt_mel := gsub("^p\\.[A-Za-z]{3}[0-9]+([A-Za-z]{3})/c\\..*","\\1", aa_change_mel)]
# filtered_dt[, aa_pos_mel := gsub(".*?([0-9]+).*", "\\1", aa_change_mel)]              # 795

filtered_dt[, aa_ref_sim := gsub("^p\\.([[:alpha:]]{3}).*", "\\1" , aa_change_sim)]   
filtered_dt[, aa_alt_sim := gsub("^p\\.[A-Za-z]{3}[0-9]+([A-Za-z]{3})/c\\..*","\\1", aa_change_sim)]

# filtered_dt[, aa_alt_sim := gsub(".*([[:alpha:]]{3})$", "\\1" , aa_change_sim)]       
# filtered_dt[, aa_pos_sim := gsub(".*?([0-9]+).*", "\\1", aa_change_sim)] 

filtered_dt[, c("aa_change_mel", "aa_change_sim") := NULL] # remove this column

##### extract codon change: ref and alt codon mel and sim : "aGc/aTc"
filtered_dt[, codon_ref_mel := gsub("^([[:alpha:]]{3})/.*", "\\1", codon_change_mel)] # aGc
filtered_dt[, codon_alt_mel := gsub(".*/([[:alpha:]]{3})$", "\\1", codon_change_mel)] # aTc

filtered_dt[, codon_ref_sim := gsub("^([[:alpha:]]{3})/.*", "\\1", codon_change_sim)] # aGc
filtered_dt[, codon_alt_sim := gsub(".*/([[:alpha:]]{3})$", "\\1", codon_change_sim)] # aTc

filtered_dt[, c("codon_change_mel", "codon_change_sim") := NULL]

##### add empty classification column
filtered_dt[, classification := NA_character_] # add empty classification colun

final_dt <- filtered_dt[, .(
  chr, pos, classification,
  variant.id_mel, variant.id_sim,
  n_samps_mel, n_samps_sim, 
  ref_mel, alt_mel, ref_sim, alt_sim, 
  codon_ref_mel, codon_alt_mel, codon_ref_sim, codon_alt_sim,
  aa_ref_mel, aa_alt_mel, aa_ref_sim, aa_alt_sim,
  af_mel, maf_mel, af_sim, maf_sim,
  effect_mel, effect_sim, 
  gene_id_fbgn, gene_mel, gene_sim  
)]

# ##### set column order 
# setcolorder(
#   filtered_dt,
#   c(
#     "chr", "pos",
#     "classification",
#     "ref_mel", "alt_mel", "ref_sim", "alt_sim",
#     "nt_change_mel", "nt_change_sim", 
#     "codon_ref_mel", "codon_alt_mel", "codon_ref_sim", "codon_alt_sim",
#     "aa_ref_mel", "aa_alt_mel", "aa_ref_sim", "aa_alt_sim",
#     "af_mel", "af_sim",
#     "effect_mel", "effect_sim",
#     "gene_mel", "gene_sim",
#     "gene_id_mel", "gene_id_sim"
#   )
# )

# check column names
names(filtered_dt) 

saveRDS(final_dt, out_rds)
# saveRDS(filtered_dt, out_rds)
message("saved clean rds to ", out_rds)

# subset_table <- filtered_dt[1:500, ]
# fwrite(subset_table, out_csv)
# message("saved first 500 rows to csv at ", out_csv)



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
 