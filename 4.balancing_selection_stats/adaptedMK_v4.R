library(SeqArray)
library(data.table)
library(doMC)
registerDoMC(16)
library(stringr)

poolseq_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/poolseq_dest/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.shared.bothMelSim.classed.fbgn.rds")
inbred_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.11_2026.shared.bothMelSim.0.05MAFfilter.classed.geneID.rds")
# gene_col <- "gene_id_fbgn_fromGFF"

mel_gff <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff"


cand_classes <- c("A", "B", "F", "G", "O", "P", "X", "Y")
tsp_classes <- c("A", "B")
conv_classes <- c("F", "G", "O", "P", "X", "Y")

NS_tsp <- c("A")
NS_conv <- c("F", "O", "X")
NS_shared <- c("A", "C", "E", "F", "H", "J", "K", "N", "O", "Q", "S", "T", "W", "X")

gff <- fread(
        cmd = paste("grep -v '^#'", shQuote(mel_gff)), 
        sep="\t", header=F, 
        col.names=c(
            "chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"
    ))

# map chr names to dm6 normal names 
chr_map <- gff[ feature == "region" & grepl("genome=chromosome", attributes),
    .(chr_gff = chr, chr = sub(".*Name=([^;]+);.*", "\\1", attributes))]

gff[chr_map, chr_dm6 := i.chr, on = .(chr = chr_gff)]
unique(gff[!is.na(chr_dm6), .(chr, chr_dm6)]) # X, 2L, 2R, 3L, 3R, 4, Y
cds_gff <- gff[feature == "CDS" & !is.na(chr_dm6), .(chr_dm6, start, end, strand, phase, attributes)]
cds_gff[, transcript_id := str_extract(attributes, "(?<=Parent=rna-)[^;,]+")]
cds_gff[, gene_name := str_extract(attributes, "(?<=;gene=)[^;]+")]
cds_gff[, gene_id_fbgn := str_extract(attributes, "FBgn[0-9]+")]


poolseq_dt[
    cds_gff,
    test_dm6_geneID_fbgn := i.gene_id_fbgn,
    on = .(
        chr = chr_dm6,
        pos >= start,
        pos <= end
    )
]



asymptotic_MKlike_stats <- function(Ps, Pns, SPs, SPns, pseudo=0, min_count=0) {
    Ps <- Ps + pseudo
    Pns <- Pns + pseudo
    SPs <- SPs + pseudo
    SPns <- SPns + pseudo

     if (Ps > min_count && Pns > min_count && SPs > min_count && SPns > min_count) {
        alpha <- 1 - (Ps * SPns) / (Pns * SPs)
        OR <- 1 - alpha 
        data.table(
            alpha = alpha, 
            OR = OR,
            logOR = log(OR)
        )
    } else {
        data.table(
            alpha = NA_real_,
            OR = NA_real_,
            logOR = NA_real_
        )
    }
}

