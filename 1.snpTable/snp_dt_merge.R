library(SeqArray)
library(data.table)
library(foreach)
library(doMC)

######################################################################

# ### Pool seq files ###

out_dir <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/"
out_csv <- paste0(out_dir, "all_quality_variants_merge_unfilt_500test.csv")
out_rds <- paste0(out_dir, "all_quality_variants_merge_unfilt.rds")
if(!file.exists(out_csv)) file.create(out_csv)
if(!file.exists(out_rds)) file.create(out_rds)

mel_snp_rds <- paste0(out_dir, "species_rdsFiles/mel_filtered_eff_snp_dt.rds")
mel_snp_dt <- readRDS(mel_snp_rds)
sim_snp_rds <- paste0(out_dir, "species_rdsFiles/sim_filtered_eff_snp_dt.rds")
sim_snp_dt <- readRDS(sim_snp_rds)

#######################################################################

### INBRED ###
out_dir <- "/scratch/ejy4bu/drosophila/inbred/snpDT/"
sim_snp_dt <- readRDS(paste0(out_dir, "dsim3.signor.snp_dt_SynMissense.rds"))
mel_snp_dt <- readRDS(paste0(out_dir, "DGRP2.source_BCM-HGSC.dm6.snp_dt_SynMissense.rds"))

shared_table <- merge(mel_snp_dt, sim_snp_dt, by = c("chr", "pos"), suffixes = c("_mel", "_sim"), all=T)
out_rds <- paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.rds")



# ### test on chromosome 2L
# # filter by 2L chromosome for a smaller subset to test merge on for csv readable
# mel_snp_dt <- mel_snp_dt[chr == "2L"]
# message(nrow(mel_snp_dt), " mel 2L variants")
# sim_snp_dt <- sim_snp_dt[chr == "2L"]
# message(nrow(sim_snp_dt), " sim 2L variants")

# NOTE: this is a union merge... keep ALL variants. no filtering yet
shared_table <- merge(mel_snp_dt, sim_snp_dt, by = c("chr", "pos"), suffixes = c("_mel", "_sim"), all=T)
message(nrow(shared_table), " total variants")


library(AnnotationDbi)
library(org.Dm.eg.db)
gene_map <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys = unique(shared_table$gene_mel),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "FLYBASE")
)

gene_map <- as.data.table(gene_map)
setnames(gene_map, c("SYMBOL", "FLYBASE"), c("gene_mel", "gene_id_fbgn"))

gene_map <- unique(gene_map[, .(gene_mel, gene_id_fbgn)])

shared_table[gene_map, gene_id_fbgn := i.gene_id_fbgn, on = "gene_mel"]

shared_table[, aa_ref_mel := gsub("^p\\.([[:alpha:]]{3}).*", "\\1" , aa_change_mel)]   # Ser
shared_table[, aa_alt_mel := gsub("^p\\.[A-Za-z]{3}[0-9]+([A-Za-z]{3})/c\\..*","\\1", aa_change_mel)]

shared_table[, aa_ref_sim := gsub("^p\\.([[:alpha:]]{3}).*", "\\1" , aa_change_sim)]   
shared_table[, aa_alt_sim := gsub("^p\\.[A-Za-z]{3}[0-9]+([A-Za-z]{3})/c\\..*","\\1", aa_change_sim)]

shared_table[, aa_pos_mel := gsub(".*?([0-9]+).*", "\\1", aa_change_mel)]
shared_table[, aa_pos_sim := gsub(".*?([0-9]+).*", "\\1", aa_change_sim)]

shared_table[, codon_ref_mel := gsub("^([[:alpha:]]{3})/.*", "\\1", codon_change_mel)] # aGc
shared_table[, codon_alt_mel := gsub(".*/([[:alpha:]]{3})$", "\\1", codon_change_mel)] # aTc

shared_table[, codon_ref_sim := gsub("^([[:alpha:]]{3})/.*", "\\1", codon_change_sim)] # aGc
shared_table[, codon_alt_sim := gsub(".*/([[:alpha:]]{3})$", "\\1", codon_change_sim)] # aTc



mel_dt <- shared_table[!is.na(ref_mel), .(
        chr, pos, 
        ref_mel, alt_mel, 
        codon_ref_mel, codon_alt_mel, 
        gene_id_fbgn, gene_mel
    )
]

sim_dt <- shared_table[!is.na(ref_sim), .(
        chr, pos, 
        ref_sim, alt_sim, 
        codon_ref_sim, codon_alt_sim,
        gene_sim
    )
]

setorder(mel_dt, chr, pos)
mel_dt[, dist_prev_mel := pos - data.table::shift(pos), by = chr]
mel_dt[, dist_next_mel := data.table::shift(pos, type = "lead") - pos, by = chr]

setorder(sim_dt, chr, pos)
sim_dt[, dist_prev_sim := pos - data.table::shift(pos), by = chr]
sim_dt[, dist_next_sim := data.table::shift(pos, type = "lead") - pos, by = chr]

shared_table <- merge(shared_table, mel_dt[, .(chr, pos, dist_prev_mel, dist_next_mel)], by=c("chr", "pos"))
shared_table <- merge(shared_table, sim_dt[, .(chr, pos, dist_prev_sim, dist_next_sim)], by=c("chr", "pos"))


unmapped_genes <- unique(
    shared_table[!is.na(variant.id_mel) & is.na(gene_id_fbgn), gene_mel]
)

cg_keys <- AnnotationDbi::keys(
    org.Dm.eg.db,
    keytype = "FLYBASECG"
)

unmapped_cg <- intersect(unmapped_genes, cg_keys)

length(unmapped_genes)
length(unmapped_cg)
cg_map <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys = unmapped_cg,
    keytype = "FLYBASECG",
    columns = c("FLYBASECG", "FLYBASE")
)

cg_map <- as.data.table(cg_map)

setnames(
    cg_map,
    c("FLYBASECG", "FLYBASE"),
    c("gene_mel", "gene_id_fbgn")
)

cg_map <- unique(cg_map[, .(gene_mel, gene_id_fbgn)])

message("saving rds to ", out_rds)
saveRDS(shared_table, out_rds)

message("complete. ", nrow(shared_table), " variants written.")

subset_table <- shared_table[1:500, ]
fwrite(subset_table, out_csv)
message("saved first 500 rows to csv at ", out_csv)
