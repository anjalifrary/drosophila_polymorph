library(data.table)
library(ggplot2)
library(foreach)
library(SeqArray)


# nlp file from Alan - contains only mel
load("/project/berglandlab/multispecies_endemism/data/collectiveAnalysis_version3/Drosophila_melanogaster.10_03_2026.nlpTable.paramask.genmap.busco.Rdata")
mel_nlp <- nlp 
load("/project/berglandlab/multispecies_endemism/data/collectiveAnalysis_version3/Drosophila_simulans.10_03_2026.nlpTable.paramask.genmap.busco.Rdata")
sim_nlp <- nlp
rm(nlp)    

variants <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed_geva.rds")

sim_nlp[, chr:=sub("^sim_", "", chr)]

sim_nlp <- merge(sim_nlp, variants[, c("chr", "pos", "codon_start_pos", "classification", "af_mel", "af_sim")], by=c("chr", "pos"), all.x=T)
# gives 1902 variants where classification is not empty ???? 
mel_nlp <- merge(mel_nlp, variants[, c("chr", "pos", "codon_start_pos", "classification", "af_mel", "af_sim")], by=c("chr", "pos"), all.x=T)

nrow(sim_nlp[!is.na(classification)])
nrow(mel_nlp[!is.na(classification)])

mel_nlp <- mel_nlp[classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y")]
nrow(mel_nlp) # 112666

sim_nlp <- sim_nlp[classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y")]
nrow(sim_nlp) # 985


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

mel_dt_5 <- merge(mel_dt_5, variants[, .(chr, pos, codon_start_pos)], by = c("chr", "pos"), all.x = TRUE)
sim_dt_5 <- merge(sim_dt_5, variants[, .(chr, pos, codon_start_pos)], by = c("chr", "pos"), all.x = TRUE)

voi <- mel_dt_5[, .(chr, pos, classification, codon_start_pos, poly_af, mel = 1L, sim = 0L)]
voi[classification %in% c("A","B","F","G","O","P"), sim := 1L]

xy_sim <- merge(
  mel_dt_5[classification %in% c("X","Y"), .(chr, codon_start_pos)],
  variants[classification %in% c("X","Y"), .(chr, pos, codon_start_pos, classification)],
  by = c("chr", "codon_start_pos")
)
xy_sim[, `:=`(mel = 0L, sim = 1L)]

xy_sim <- xy_sim[, .(chr, pos, classification, mel, sim)]

voi_final <- rbindlist(list(voi[, .(chr, pos, classification, mel, sim)], xy_sim[, .(chr, pos, classification, mel, sim)]),
  use.names = TRUE
)
table(voi_final$classification)

saveRDS(voi_final, "/scratch/ejy4bu/drosophila/gowinda/maf_filter_mel5/voi_abfgopxy.rds")


# merge tsp file on chr and pos, include classification column
# nlp <- merge(nlp, tsp[,c("chr", "pos", "classification")], by=c("chr", "pos"), all.x=T)
dt <- merge(tsp, nlp[, c("chr", "pos", "busco", "genmap_score", "nLocales_poly", "global_af", "poly_af", "poly_maf")], by=c("chr", "pos"), all.x=T)

age <- fread("/scratch/ejy4bu/drosophila/AlleleAges.VA.cm_GEVA.txt")

# merge nlp with age on chr and pos 
age[,chr:=tstrsplit(id, "\\.")[[1]]]
age[,pos:=position]
dt <- merge(dt, age, by=c("chr", "pos"), all.x=T)

all_var <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/merged_tables/quality/all_quality_variants_merge_unfilt.rds")
# column names: chr, pos, variant.id
dt <- merge(dt, all_var[, c("chr", "pos", "variant.id_mel", "variant.id_sim")], by=c("chr", "pos"), all.x=T)

mel_meta <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/dest_v2.samps_24Aug2024.xa.csv")
setnames(mel_meta, "sampleId", "sample.id")

# column names: sampleId, lat, long
sim_meta <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/simulans_pooled.meta.use.csv")
setnames(sim_meta, "identifier", "sample.id")

mel_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.eff.gds")
sim_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.eff.gds")

mel_samples <- data.table( sample.id = seqGetData(mel_gds, "sample.id"))
sim_samples <- data.table( sample.id = seqGetData(sim_gds, "sample.id"))

mel_samples <- merge(mel_samples, mel_meta, by = "sample.id", all.x = TRUE)
mel_samples[, sample_index := .I]

tsp_ids <- dt[, variant.id_mel]
seqSetFilter(mel_gds, variant.id = tsp_ids)
geno <- seqGetData(mel_gds, "genotype")


# now what??? 
# what about sim ?? 
# how to visualize the geographic mapping?

# adding MAF filter 