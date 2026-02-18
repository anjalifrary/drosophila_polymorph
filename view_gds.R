# .libPaths(c("~/Rlibs", .libPaths()))
library(SeqArray)
library(SeqVarTools)

out_file <- "/scratch/ejy4bu/drosophila/gds_analysis/view_gds.txt"
if(!file.exists(out_file)) file.create(out_file)

mel_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.gds"
sim_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds"

mel_gds <- seqOpen(mel_file)
sim_gds <- seqOpen(sim_file)

# view stats
mel_gds
sim_gds 

# reference allele freq of each variant
mel_freq <- seqAlleleFreq(mel_gds)
sim_freq <- seqAlleleFreq(sim_gds)
head(mel_freq)
summary(mel_freq)
# hist(mel_freq, breaks=50)
mel_df <- data.frame(freq = mel_freq)
sim_df <- data.frame(freq = sim_freq)

### histogram of allele frequency

## allele freq = (# alt alleles observed)/(# total alleles observed)
# near 1 => fixed alleles
# near 0.5 => SNP
# near 0 => rare allele

library(ggplot2)

# mel:
p <- ggplot(mel_df, aes(x=freq)) + 
    geom_histogram(binwidth=0.05) +
    labs(title = "Distribution of Mel Allele Frequencies",
    x = "Allele Frequency",
    y = "Count")
ggsave("/scratch/ejy4bu/drosophila/gds_analysis/mel_freq_histogram.jpg", plot = p, width=8,height=6)

# sim
p <- ggplot(sim_df, aes(x=freq)) + 
    geom_histogram(binwidth=0.05) +
    labs(title = "Distribution of Sim Allele Frequencies",
    x = "Allele Frequency",
    y = "Count")
ggsave("/scratch/ejy4bu/drosophila/gds_analysis/sim_freq_histogram.jpg", plot = p, width=8,height=6)
