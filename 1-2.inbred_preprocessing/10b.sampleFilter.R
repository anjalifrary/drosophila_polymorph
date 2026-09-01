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
library(data.table)


outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"
## mel
mel_gds <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.gds"
mel_snp_gds <- paste0(outdir, "DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.snpRelate.gds")
## sim
sim_gds <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.gds"
sim_snp_gds <- paste0(outdir, "dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.snpRelate.gds")

# meta_file <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/DGRP2.source_BCM-HGSC.dm6.csv"
meta_file <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/signor.dsim3.sampleFilt.csv"

# array_genofile <- seqOpen(mel_gds)
# relate_genofile <- snpgdsOpen(mel_snp_gds)

array_genofile <- seqOpen(sim_gds)
# relate_genofile <- snpgdsOpen(sim_snp_gds)

##################################################################
### MISSINGNESS
##################################################################


# calculate missing rate using SeqArray:
sample_missing_array <- seqMissing(array_genofile, per.variant=FALSE)
samples <- seqGetData(array_genofile, "sample.id")
dt <- data.frame(sample.id=samples, missing.seqarray=sample_missing_array)

metadata <- as.data.table(dt)
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



##################################################################
### AVG RD
##################################################################


######## MEL #########

# for mel, RD is stored as OR + SR 
# $format (annotation/format/OR)
# 3 OR      1 Integer   number of opposing reads
# 4 SR      1 Integer number of supporting reads

rd_sum <- numeric(length(samples))
rd_count <- numeric(length(samples))

seqApply(array_genofile, c("annotation/format/OR", "annotation/format/SR"),
    function(x) {
        or <- x[[1]]
        sr <- x[[2]]
        rd <- or + sr
        rd_sum <<- rd_sum + rowSums(rd, na.rm = TRUE) # sum of rd across all rows
        rd_count <<- rd_count + rowSums(!is.na(rd)) # number of rows that not NA for rd 
        NULL
    },
    margin = "by.variant")

avg_rd <- rd_sum / rd_count

rd_table <- data.frame(sample.id = samples, avg.RD = avg_rd)
summary(rd_table$avg.RD)
head(rd_table)

hist(rd_table$avg.RD, breaks = 50)

metadata <- merge(metadata, rd_table, by="sample.id")

######## SIM #########
# sim gds does not have OR or SR. instead uses DP 
# DP, 1, Integer, Approximate read depth (reads with MQ=255 or with bad mates are filtered)

rd_sum <- numeric(length(samples))
rd_count <- numeric(length(samples))

seqApply(array_genofile, "annotation/format/DP",
    function(dp) {
        rd_sum <<- rd_sum + rowSums(dp, na.rm = TRUE)
        rd_count <<- rd_count + rowSums(!is.na(dp))
        NULL
    },
    margin = "by.variant"
)

avg_rd <- rd_sum / rd_count

rd_table <- data.frame(sample.id = samples, avg.RD = avg_rd)
summary(rd_table$avg.RD)
head(rd_table)

hist(rd_table$avg.RD, breaks = 50)

metadata <- merge(metadata, rd_table, by="sample.id")



##################################################################
### Proportion heterozygous snps 
##################################################################

# num heterozygous snps / num non-missing snps
# samples <- seqGetData(array_genofile, "sample.id")

het_count <- numeric(length(samples))
called_count <- numeric(length(samples))

# seqGetData(array_genofile, "genotype") gives 3D array (2 alleles, 205 samples, xx variants)

seqApply(array_genofile, "genotype",
    function(x) {
        het <- x[1, ] != x[2, ] # all samples x all variants, comparing alleles 1 and 2 and selecting those not equal
        called <- !is.na(x[1, ]) & !is.na(x[2, ]) # all sampls x all variants, only counting nonmissing
        het_count <<- het_count + ifelse(het & called, 1, 0) # update het_count outside of function
        called_count <<- called_count + ifelse(called, 1, 0) # update called_count 
        NULL
    },
    margin = "by.variant" # process genotype in chunks by variant instead of simultaneously 
)
prop_heterozygous <- het_count / called_count

het_table <- data.frame(
    sample.id = samples,
    het.count = het_count,
    called.count = called_count,
    prop.heterozygous = prop_heterozygous
)

head(het_table)
summary(het_table$prop.heterozygous)

# merge het_table on metadata
metadata <- merge(metadata, het_table, by="sample.id")



fwrite(metadata, meta_file)
seqClose(array_genofile)
# snpgdsClose(relate_genofile)