library(data.table)

inbred <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.candidatesABFGOPXY.classed.MAF.rds")
poolseq <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/poolseq_dest/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.candidatesABFGOPXY.classed.maf.nlp.xtx.geva.rds")


# unique genomic positions in each dataset
inbred_pos <- unique(inbred[, .(chr_dm6, pos_dm6)])
poolseq_pos <- unique(poolseq[, .(chr_dm6, pos_dm6)])

# intersection
shared_pos <- merge(
    inbred_pos,
    poolseq_pos,
    by = c("chr_dm6", "pos_dm6")
)
nrow(shared_pos)

shared <- merge(
    inbred[, .(chr_dm6, pos_dm6, classification)],
    poolseq[, .(chr_dm6, pos_dm6, classification)],
    by = c("chr_dm6", "pos_dm6"), suffixes=c("_inbred", "_poolseq"),
    all = FALSE
)
table(shared$type, useNA = "ifany")

no_match <- shared[classification_inbred != classification_poolseq]
no_match[, group_inbred := ifelse(no_match$classification_inbred%in%c("A", "B"), "tsp", "conv")]
no_match[, group_poolseq := ifelse(no_match$classification_poolseq%in%c("A", "B"), "tsp", "conv")]