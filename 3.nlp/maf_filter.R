library(data.table)
library(ggplot2)
library(foreach)
library(SeqArray)

# sim file
load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_simulans.17_06_2026.nlpTable.Rdata")
# load("/project/berglandlab/multispecies_endemism/data/collectiveAnalysis_version3/Drosophila_simulans.10_03_2026.nlpTable.paramask.genmap.busco.Rdata")
### this sim file is in simulans genospace, wrong coordinates
sim_nlp <- nlp 
rm(nlp)

# mel file
# load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_melanogaster.17_06_2026.nlpTable.Rdata")
load("/project/berglandlab/multispecies_endemism/data/collectiveAnalysis_version3/Drosophila_melanogaster.10_03_2026.nlpTable.paramask.genmap.busco.Rdata")
# confirm these nlp files are the same... (excluding extra columns)
mel_nlp <- nlp
rm(nlp)   

cand_classes <- c("A", "B", "F", "G", "O", "P", "X", "Y")

# all variant table, classed, no maf filter yet:
variants <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed_geva.rds")

sim_nlp <- merge(sim_nlp, variants[, c("chr", "pos", "codon_start_pos", "classification")], by=c("chr", "pos"), all.x=T)
mel_nlp <- merge(mel_nlp, variants[, c("chr", "pos", "codon_start_pos", "classification")], by=c("chr", "pos"), all.x=T)

nrow(sim_nlp[!is.na(classification)])
nrow(mel_nlp[!is.na(classification)])

nrow(mel_nlp) 
nrow(sim_nlp) 

nrow(mel_nlp[classification%in%(cand_classes)]) # 112666
nrow(sim_nlp[classification%in%(cand_classes)]) # 985

# filter for polymorphisms in which nlp > 10 *** check with Alan
sim_nlp <- sim_nlp[nLocales_poly>10]
mel_nlp <- mel_nlp[nLocales_poly>10]

# recheck metrics
nrow(mel_nlp) 
nrow(sim_nlp) 

nrow(mel_nlp[classification%in%(cand_classes)]) 
nrow(sim_nlp[classification%in%(cand_classes)]) 



# what allele frequency to be considered a polymorphism ?? 
nrow(mel_nlp[poly_af>0.05 & poly_af<0.95]) # 56397
nrow(sim_nlp[poly_af>0.05 & poly_af<0.95]) # 821

nrow(mel_nlp[poly_af>0.03 & poly_af<0.97]) # 90185
nrow(sim_nlp[poly_af>0.03 & poly_af<0.97]) # 931

nrow(mel_nlp[poly_af>0.01 & poly_af<0.99]) # 112178
nrow(sim_nlp[poly_af>0.01 & poly_af<0.99]) # 985


# 5% MAF filter:
mel_dt_5 <- mel_nlp[poly_af>0.05 & poly_af<0.95]
sim_dt_5 <- sim_nlp[poly_af>0.05 & poly_af<0.95]

mel_dt_5 <- merge(mel_dt_5[!is.na(ref_mel), .(chr, pos, variant, nLocales_poly, global_af, poly_af, poly_maf)], variants, by=c("chr", "pos"), all.x=T)
sim_dt_5 <- merge(sim_dt_5[!is.na(ref_sim), .(chr, pos, variant, nLocales_poly, global_af, poly_af, poly_maf)], variants, by=c("chr", "pos"), all.x=T)


voi <- rbindlist(list(
    mel_dt_5, 
    sim_dt_5
))

names(voi)


voi_candidates <- voi[classification%in%(cand_classes)]

saveRDS(voi, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed_geva_MAF5.rds")
saveRDS(voi_candidates, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed_geva_MAF5_ABFGOPXY.rds")



##### trash #####
# from before I made the proper simulans nlp file

# ### to get only variants required for gowinda analysis:
# mel_dt_5 <- merge(mel_dt_5, variants[, .(chr, pos, codon_start_pos)], by = c("chr", "pos"), all.x = TRUE)
# sim_dt_5 <- merge(sim_dt_5, variants[, .(chr, pos, codon_start_pos)], by = c("chr", "pos"), all.x = TRUE)

# voi <- mel_dt_5[, .(chr, pos, classification, codon_start_pos, poly_af, mel = 1L, sim = 0L)]
# # voi[classification %in% c("A","B","F","G","O","P"), sim := 1L]

# xy_sim <- merge(
#   mel_dt_5[classification %in% c("X","Y"), .(chr, codon_start_pos)],
#   variants[classification %in% c("X","Y"), .(chr, pos, codon_start_pos, classification)],
#   by = c("chr", "codon_start_pos")
# )
# xy_sim[, `:=`(mel = 0L, sim = 1L)]

# xy_sim <- xy_sim[, .(chr, pos, classification, mel, sim)]

# voi_final <- rbindlist(list(voi[, .(chr, pos, classification, mel, sim)], xy_sim[, .(chr, pos, classification, mel, sim)]),
#   use.names = TRUE
# )
# table(voi_final$classification)

# saveRDS(voi_final, "/scratch/ejy4bu/drosophila/gowinda/maf_filter_mel5/voi_abfgopxy.rds")

# ### to get variants of interest file with MAF 5% in mel filter:
# mel_dt_5 <- merge(mel_dt_5[, .(chr, pos, variant, nLocales_poly, global_af, poly_af, poly_maf)], variants, by=c("chr", "pos"), all.x=T)
# sim_dt_5 <- merge(sim_dt_5[, .(chr, pos, variant, nLocales_poly, global_af, poly_af, poly_maf)], variants, by=c("chr", "pos"), all.x=T)

# voi <- rbindlist(list(
#     mel_dt_5[classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y")], 
#     sim_dt_5[classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y")]
# ), use.names=T)

# voi_candidates <- voi[classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y")]

# saveRDS(voi, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed_geva_melMAF5.rds")
