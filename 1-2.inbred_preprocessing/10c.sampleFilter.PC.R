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

library(seqArray)
library(SNPRelate)


# snpgdsPCA
gds <- "/"

genofile <- seqOpen(gds)
samps <- seqGetData(genofile, var.name = "sample.id")


# get sample metadata, make new column for this 


# pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

# connor's version:
# ccm_pca <- snpgdsPCA(genofile, 
#                      autosome.only = F, 
#                      sample.id = as.character(fin$Sample),
#                      num.thread = 10, 
#                      snp.id = unique(snps$variant.id),
#                      maf=0.01)

saveRDS(pca, ".rds")


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
