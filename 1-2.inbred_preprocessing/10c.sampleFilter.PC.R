### tbd 
# if removing individuals, merge with 10b.sampleFilter.R

# using package snp relate 
# to visualize population structure 

# principle components - correlation
# PC1 , PC2, plot against one another
# expecting 1 "blob" - can probably remove outlier individuals which may be low coverage etc 

# connor:
# https://github.com/connor122721/SharedPolymorphismsDaphnia/blob/864eda3b218e1bd2880907a651d791165021261e/Figures/Figure_1_Phylogeny/Scripts/Run_Multi_Species_PCA.R#L35

# manual
# https://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html#principal-component-analysis-pca

library(SeqArray)
library(SNPRelate)


# snpgdsPCA
outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"
## mel
mel_gds <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.gds"
mel_snp_gds <- paste0(outdir, "DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.snpRelate.gds")
## sim
sim_gds <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.gds"
sim_snp_gds <- paste0(outdir, "dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.snpRelate.gds")


# get SNPRelate gds files, (convert from SeqArray gds file formats)
# seqGDS2SNP(mel_gds, mel_snp_gds)
genofile <- snpgdsOpen(mel_snp_gds)

# seqGDS2SNP(sim_gds, sim_snp_gds)
# genofile <- snpgdsOpen(sim_snp_gds)

genofile

# genofile <- seqOpen(mel_gds)
# genofile <- seqOpen(sim_gds)
# samps <- seqGetData(genofile, var.name = "sample.id")

samples <- read.gdsn(index.gdsn(genofile, "sample.id"))

head(samples)

length(samples) # should be 205 for mel; 183 for sim

# LD pruning?? to get relatively independent snp set for PCA?
# snpgdsLDpruning()
set.seed(1000)
snpset <- snpgdsLDpruning(
    genofile,
    # maf=0.01,
    # missing.rate = 0.01,
    method="corr",
    autosome.only = FALSE,
    slide.max.bp = 500000,
    ld.threshold=0.2,
    num.thread = 10
)
# optional params:
# missing.rate? to use the SNPs with "<= missing.rate" only; if NaN, no missing threshold
# maf? to use the SNPs with ">= maf" only; if NaN, no MAF threshold
# slide.max.bp? the maximum basepairs in the sliding window

str(snpset)
names(snpset)
snpIDs <- unlist(unname(snpset))

length(snpIDs) # orig 2.95M snps

sapply(snpset, length)

pca <- snpgdsPCA(
    genofile,
    # snp.id = snpIDs,
    maf = 0.01,
    # missing.rate = 0.01,
    autosome.only = FALSE,
    num.thread = 10
)

pc.percent <- pca$varprop * 100
head(round(pc.percent, 2))

table <- data.frame(
    sample.id = pca$sample.id,
    PC1 = pca$eigenvect[, 1],
    PC2 = pca$eigenvect[, 2]
)

head(table)

library(ggplot2)

ggplot(table, aes(x = PC1, y = PC2)) +
    geom_point(size = 1) +
    theme_classic() +
    labs(
        x = paste0("PC1 (", round(pc.percent[1], 2), "%)"),
        y = paste0("PC2 (", round(pc.percent[2], 2), "%)"),
        title = "Mel PCA"
    )

# get missingness per sample: 
sample_missing <- snpgdsSampMissRate(genofile)
table$missing.rate <- sample_missing

# color by missingness
ggplot(table, aes(PC1, PC2, color = missing.rate)) +
    geom_point(size = 1) +
    theme_classic() + 
    labs(
        x = paste0("PC1 (", round(pc.percent[1], 2), "%)"),
        y = paste0("PC2 (", round(pc.percent[2], 2), "%)"),
        title = "Mel PCA - unpruned & color-coded",
        color = "Missing rate"
    )

## identify outliers 
library(plotly)

p <- ggplot(table, aes(PC1, PC2, text = sample.id, color = missing.rate)) +
    geom_point(size = 1) +
    theme_classic()

ggplotly(p, tooltip = "text")

# plot(
#     table$PC1,
#     table$PC2,
#     xlab = paste0("PC1 (", round(pc.percent[1], 2), "%)"),
#     ylab = paste0("PC2 (", round(pc.percent[2], 2), "%)")
# )

# outlier_ids <- identify(
#     table$PC1,
#     table$PC2,
#     labels = table$sample.id
# )

# table$sample.id[outlier_ids]



snpgdsClose(genofile)

# saveRDS(pca, paste0(outdir, "/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.PCA.rds"))


# get gds in format for SNPRelate
# get pca1 and pca2
# find correlation + outliers + filter

# get sample metadata, make new column for this 


# pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

# connor's version:
# ccm_pca <- snpgdsPCA(genofile, 
#                      autosome.only = F, 
#                      sample.id = as.character(fin$Sample),
#                      num.thread = 10, 
#                      snp.id = unique(snps$variant.id),
#                      maf=0.01)

# saveRDS(pca, ".rds")


# The code below shows how to calculate the percent of variation is accounted for by the top principal components. It is clear to see the first two eigenvectors hold the largest percentage of variance among the population, although the total variance accounted for is still less the one-quarter of the total.

# # variance proportion (%)
# pc.percent <- pca$varprop*100
# head(round(pc.percent, 2))

# # In the case of no prior population information,

# # make a data.frame
# tab <- data.frame(sample.id = pca$sample.id,
#     EV1 = pca$eigenvect[,1],    # the first eigenvector
#     EV2 = pca$eigenvect[,2],    # the second eigenvector
#     stringsAsFactors = FALSE)
# head(tab)

# plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
