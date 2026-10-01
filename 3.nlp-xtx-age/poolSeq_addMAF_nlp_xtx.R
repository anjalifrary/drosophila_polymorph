library(data.table)
library(ggplot2)
library(foreach)
library(SeqArray)

shared_dt <- readRDS("/scratch/ejy4bu/drosophila/DEST_remake/snpDT/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.shared.bothMelSim.classed.rds")

# shared_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/poolseq_dest/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.shared.bothMelSim.classed.rds")


final_dt_save <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.shared.bothMelSim.classed.maf.nlp.xtx.geva.rds"
sharedOnly_dt <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.shared.classed.maf.nlp.xtx.geva.rds"
candidates_dt <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.candidatesABFGOPXY.classed.maf.nlp.xtx.geva.rds"


# xtx
load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/xtx_c2.Rdata")
xtx <- xc
rm(xc)

# nlp, poly_AF, global_AF
load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_melanogaster.11_08_2026.nlpTable.paramask.genmap.busco.repeatmask.wmdust.Rdata")
mel_nlp <- nlp  
rm(nlp)
load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_simulans.17_06_2026.nlpTable.Rdata")
sim_nlp <- nlp
rm(nlp)


age <- fread("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/AlleleAges.VA.cm_GEVA.txt")
age[,chr:=tstrsplit(id, "\\.")[[1]]]
age[,pos:=position]

# merge with xtx
shared_dt <- merge(shared_dt, xtx[, .(chr, pos, XtXst)], by.y=c("chr", "pos"), by.x=c("chr_dm6", "pos_dm6"), all.x=T)

# merge with nlp 
shared_dt <- merge(shared_dt, mel_nlp[, .(chr, pos, nLocales_poly, global_af, poly_af, busco)],
    by.x=c("chr_dm6", "pos_dm6"), by.y=c("chr", "pos"), all.x=T)
setnames(shared_dt, c("nLocales_poly", "global_af", "poly_af", "busco"), 
    c("nLocales_poly_mel", "global_af_mel", "poly_af_mel", "busco_mel"))

shared_dt <- merge(shared_dt, sim_nlp[, .(chr, pos, nLocales_poly, global_af, poly_af)],
    by.x=c("chr_dm6", "pos_dm6"), by.y=c("chr", "pos"), all.x=T)
setnames(shared_dt, c("nLocales_poly", "global_af", "poly_af"), 
    c("nLocales_poly_sim", "global_af_sim", "poly_af_sim"))


# merge with age:

shared_dt <- merge(shared_dt, age[, .(chr, pos, PostMode, PostMean, PostMedian)], by.y=c("chr", "pos"), by.x=c("chr_dm6", "pos_dm6"), all.x=T)

final_dt <- shared_dt
sharedOnly <- final_dt[!is.na(classification)]
candidates <- final_dt[classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y")]


saveRDS(final_dt, final_dt_save)

saveRDS(sharedOnly, sharedOnly_dt)

saveRDS(candidates, candidates_dt)