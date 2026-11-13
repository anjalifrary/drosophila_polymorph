library(data.table)
inbred <- readRDS("/scratch/ejy4bu/drosophila/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.shared.bothMelSim.0.05MAFfilter.classed.rds")
# poolseq <- 

snp.dt <- inbred

# calculated by adding avg-lineage-specific divergence via Langley, 2012, Genetics, where Dmel = avg(0.0384, 0.03017); Dsim = 0.02432
dMEL <- mean(c(0.0384, 0.03017)) 
dSIM <- 0.02432
dXY <- mean(c(dMEL, dSIM)) # between species Dxy

# syn polymorphism within species; calculated by Nolte, 2008, Genetics, taking averages 
# note: expected psim > pmel 
pSIM <- mean(c(0.0107, 0.0131)) # within sim 
pMEL <- mean(c(0.0114, 0.0146)) # within mel

e1 <- exp(-((dXY-max(pSIM, pMEL))/(pSIM)))
e2 <- exp(-((dXY-max(pSIM, pMEL))/(pMEL)))
prop.exp.tsp <- e1*e2
prop.exp.tsp

num.TSPcandidates <- nrow(snp.dt[classification%in%c("A", "B")])

# what should this be? total number of shared snps? like anything classed same-site?

# all same-site shared snps (includes TSP candidates, conv, divergent)
num.snps <- nrow(snp.dt[classification%in%c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")])
# num.snps <- nrow(snp.dt[!is.na(classification)]) # all shared homologous snps
# num.snps <- nrow(snp.dt) # all unique SNPs (including species-specific snps)
prop.obs.tsp <- num.TSPcandidates/num.snps

# Expected # TSPs
tsps.exp <- prop.tsp * num.snps

enrichment <- log2(num.TSPcandidates/(tsps.exp))

prop.exp.tsp # 0.07251249
tsps.exp # 782.9899
enrichment # 3.346609

prop.obs.tsp # 0.7376366
num.TSPcandidates # 7965

