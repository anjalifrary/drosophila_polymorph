############  Adding codon identity   ############################################
# don't do this anymore because I updated the datatables at the start (gds_to_sp_dt.R)

library(SeqArray)
library(data.table)

# load gds file
mel_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.eff.gds"
gds_file <- seqOpen(mel_file)
# sim_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.eff.gds"
# gds_file <- seqOpen(sim_file)

# load rds files
datatable <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/mel_snp_dt.rds")
# datatable <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/sim_snp_dt.rds")
rds_file <- readRDS(datatable)
message(nrow(rds_file), " total variants in table")


# get ids of variants already in rds file
seqSetFilter(gds_file, variant.id = rds_file$variant.id)
eff_all <- seqGetData(gds_file, "annotation/info/ANN")
annotated_ids <- seqGetData(gds_file, "variant.id")

# ann_all <- seqGetData(gds_file, "annotation/info/ANN")
# annotated_ids <- seqGetData(gds_file, "variant.id")

eff_dt <- data.table(variant.id = rep(annotated_ids, times = eff_all$length), eff = eff_all$data)
eff_split <- tstrsplit(eff_dt$eff, "\\|")
eff_dt[, codon_change := eff_split[[3]]]    # codon change. format : aAt/aCt

# keep first eff per variant
eff_dt_first <- eff_dt[, .SD[1], by = variant.id]

# merge
rds_file <- merge(rds_file, eff_dt_first[, .(variant.id, codon_change)], by="variant.id", all.x=T)

# ann_dt <- data.table(variant.id = rep(annotated_ids, times = ann_all$length), ann = ann_all$data)
# ann_dt[, effect_order := seq_len(.N), by = variant.id]
# ann_dt <- ann_dt[effect_order == 1] # confirm variants kept are 1st effect

# ann_split <- tstrsplit(ann_dt$ann, "\\|")
# ann_dt[, aa_codon := ann_split[[12]]]       # codon that codes for amino acid 
# codon_split <- tstrsplit(ann_dt$codon_change, "/")  # parse codon "cCg/cAg"
# ann_dt[, codon_ref := codon_split[[1]]]
# ann_dt[, codon_alt := codon_split[[2]]]
# ann_dt[, ann := NULL]

# rds_file <- merge(rds_file, ann_dt[, .(variant.id, codon_ref, codon_alt)], by = "variant.id", all.x=T)
saveRDS(rds_file, datatable)
message("updated: ", nrow(rds_file), " variants")
