###
## sample-level filtering:
## on avg RD, missingness, and proportion heterozygous

# use this for missingness per-snp via SNPRelate: 
# https://rdrr.io/bioc/SNPRelate/man/snpgdsSNPRateFreq.html

### to run after vcf2gds (snpEff -> liftover -> vcf2gds -> sample-level filtering x2)

# take a random snp subset  (like 100k or something)
# then per sample, extract avg RD, missing data, proportion heterozygous 
    # inbred should have low heterozygosity 


library(SeqArray)
library(SNPRelate)


outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"
## mel
mel_gds <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.gds"
mel_snp_gds <- paste0(outdir, "DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.snpRelate.gds")
## sim
sim_gds <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.gds"
sim_snp_gds <- paste0(outdir, "dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.snpRelate.gds")

array_genofile <- seqOpen(mel_gds)
relate_genofile <- snpgdsOpen(mel_snp_gds)

# array_genofile <- seqOpen(sim_gds)
# relate_genofile <- snpgdsOpen(sim_snp_gds)

##################################################################
### MISSINGNESS
##################################################################


# calculate missing rate using SeqArray:
sample_missing_array <- seqMissing(array_genofile, per.variant=FALSE)
samples <- seqGetData(array_genofile, "sample.id")
dt <- data.frame(sample.id=samples, missing.seqarray=sample_missing_array)

# # calculate missing rate using SNPRelate:
# sample_missing_snprelate <- snpgdsSampMissRate(relate_genofile, with.id = TRUE)
# # merge snprelate samples to compare
# dt$missing.snprelate <- sample_missing_snprelate[dt$sample.id]

# # look at the difference
# dt$diff <- dt$missing.seqarray - dt$missing.snprelate
# dt[order(-dt$missing.seqarray), ]

# head(dt)
# summary(dt$diff)
# max(abs(dt$diff)) # is exactly 0 - perfect.

### View missingness:
hist(dt$missing.seqarray, breaks=100)

seqClose(array_genofile)
snpgdsClose(relate_genofile)


##################################################################
### AVG RD
##################################################################
