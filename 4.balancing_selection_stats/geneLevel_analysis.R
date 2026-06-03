

library(data.table)
library(ggplot2)
library(foreach)
library(SeqArray)


# nlp file from Alan - contains only mel
load("/project/berglandlab/multispecies_endemism/data/collectiveAnalysis_version3/Drosophila_melanogaster.10_03_2026.nlpTable.paramask.genmap.busco.Rdata")

tsp <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed_geva.rds")

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

