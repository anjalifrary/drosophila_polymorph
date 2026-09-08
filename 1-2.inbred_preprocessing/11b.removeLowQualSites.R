library(SeqArray)
library(data.table)

### mel
# outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"
# mel_gds <- paste0(outdir, "DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.goodSamps.gds")
# out_gds <- paste0(outdir, "DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.goodSamps.goodSites.gds")
# meta_file <- as.data.table(read.csv("/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/DGRP2.source_BCM-HGSC.dm6.csv"))
# genofile <- seqOpen(mel_gds)
# site_rd <- readRDS(paste0(outdir, "mel_site_RD.rds"))

### sim
outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"
sim_gds <- paste0(outdir, "dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.goodSamps.gds")
out_gds <- paste0(outdir, "dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.goodSamps.goodSites.gds")
meta_file <- as.data.table(read.csv("/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/signor.dsim3.sampleFilt.csv"))
genofile <- seqOpen(sim_gds)
site_rd <- readRDS(paste0(outdir, "sim_site_RD.rds"))


#mel:
# avg_RD_min <- 5
# avg_RD_max <- 55

# keep <- site_rd[
#     avg.RD >= avg_RD_min &
#     avg.RD <= avg_RD_max
# ]

# sim:
avg_RD_min <- 3
avg_RD_max <- 15

sum_RD_min <- 500
sum_RD_max <- 2500

keep <- site_rd[
    avg.RD >= avg_RD_min &
    avg.RD <= avg_RD_max &
    sum.RD >= sum_RD_min & 
    sum.RD <= sum_RD_max
]

seqSetFilter(
    genofile,
    variant.id = keep$variant.id
)

# Check the filter
seqSummary(genofile)


seqExport(genofile,out_gds)
seqClose(genofile)



### get lower and upper bounds
Q1 <- quantile(site_rd$avg.RD, 0.25, na.rm = TRUE)
Q3 <- quantile(site_rd$avg.RD, 0.75, na.rm = TRUE)
IQR_RD <- Q3 - Q1

# 1.5 IQR was too strict maybe... 

lower <- Q1 - 2 * IQR_RD
upper <- Q3 + 2 * IQR_RD

c(lower = lower, upper = upper)

lower
upper

# filtered_dt <- site_rd[avg.RD < lower | avg.RD > upper]


### get z scores
site_rd[, avg.RD.z := 
    (avg.RD - median(avg.RD, na.rm = TRUE)) /
    mad(avg.RD, na.rm = TRUE)
]

site_rd[, sum.RD.z := 
    (sum.RD - median(sum.RD, na.rm = TRUE)) /
    mad(sum.RD, na.rm = TRUE)
]

summary(site_rd$avg.RD.z)
summary(site_rd$sum.RD.z)


### get quantiles 
quantile(
    site_rd$avg.RD,
    probs = c(0, .001, .005, .01, .02, .05, .95, .98, .99, .995, .999, 1),
    na.rm = TRUE
)

quantile(
    site_rd$sum.RD,
    probs = c(0, .001, .005, .01, .02, .05, .95, .98, .99, .995, .999, 1),
    na.rm = TRUE
)
