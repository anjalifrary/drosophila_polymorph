# calculate site-level avg RD and sum RD across all samples

args <- commandArgs(trailingOnly = TRUE)

start <- as.integer(args[1])
end <- as.integer(args[2])


library(SeqArray)
library(SNPRelate)
library(data.table)


outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"
## mel
mel_gds <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.gds"
## sim
sim_gds <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.gds"

meta_file <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/DGRP2.source_BCM-HGSC.dm6.csv"
array_genofile <- seqOpen(mel_gds)

# meta_file <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/signor.dsim3.sampleFilt.csv"
# array_genofile <- seqOpen(sim_gds)


seqResetFilter(array_genofile)

seqSetFilter(
    array_genofile,
    variant.sel = start:end,
    verbose = FALSE
)

variants <- seqGetData(array_genofile, "variant.id")
nvar <- length(variants)
print(paste0(nvar, " variants"))

### sim:
# dp <- seqGetData(array_genofile, "annotation/format/DP")

# sum.RD <- colSums(dp, na.rm = TRUE)
# called.samples <- colSums(!is.na(dp))
# avg.RD <- sum.RD / called.samples

# site_rd <- data.frame(
#     variant.id = variants,
#     sum.RD = sum.RD,
#     avg.RD = avg.RD,
#     called.samples = called.samples
# )


### mel:
OR <- seqGetData(array_genofile, "annotation/format/OR")
SR <- seqGetData(array_genofile, "annotation/format/SR")
sum.RD <- colSums(OR + SR, na.rm = TRUE)
called.samples <- colSums(!is.na(OR + SR))
# avg.RD <- colMeans(OR + SR, na.rm = TRUE)
avg.RD <- sum.RD / called.samples

avg.OR <- colMeans(OR, na.rm = TRUE)
avg.SR <- colMeans(SR, na.rm = TRUE)

site_rd <- data.frame(
    variant.id = variants,
    avg.OR = avg.OR,
    avg.SR = avg.SR,
    avg.RD = avg.RD,
    sum.RD = sum.RD,
    called.samples = called.samples
)

head(site_rd)
summary(site_rd$avg.RD)
summary(site_rd$sum.RD)

outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"
outfile <- file.path(
    outdir,
    paste0("mel_site_RD_", start, "_", end, ".rds")
)

saveRDS(site_rd, outfile)
seqClose(array_genofile)


####################################
# # combine intermediate files
# library(data.table)

# outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"
# files <- list.files(
#     outdir,
#     pattern = "^mel_site_RD_[0-9]+_[0-9]+\\.rds$",
#     full.names = TRUE
# )

# files <- sort(files)
# print(files)
# site_rd_list <- lapply(files, readRDS)

# site_rd <- rbindlist(site_rd_list)

# saveRDS(site_rd, file.path(outdir, "sim_site_RD.rds"))

# print(site_rd)
# print(dim(site_rd))