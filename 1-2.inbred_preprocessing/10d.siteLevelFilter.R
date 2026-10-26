# calculate site-level avg RD and sum RD across all samples

args <- commandArgs(trailingOnly = TRUE)

start <- as.integer(args[1])
end <- as.integer(args[2])
# start <- 1
# end <- 1000


library(SeqArray)
library(SNPRelate)
library(data.table)


outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"

### Mel - DGRP2
mel_gds <- "/scratch/ejy4bu/drosophila/inbred/gds_files/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.gds"
meta_file <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/DGRP2.source_BCM-HGSC.dm6.10_2026.csv"
outfile <- file.path(outdir, paste0("mel_site_RD_", start, "_", end, ".rds"))
array_genofile <- seqOpen(mel_gds)

### Sim - Signor
# sim_gds <- "/scratch/ejy4bu/drosophila/inbred/gds_files/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.dm6.sorted.m2.M2.gds"
# meta_file <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/signor.dsim3.sampleFilt.dm6.10_2026.csv"
# outfile <- file.path(outdir, paste0("sim_site_RD_", start, "_", end, ".rds"))
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

# ### sim:
# dp <- seqGetData(array_genofile, "annotation/format/DP")

# sum.RD <- colSums(dp, na.rm = TRUE)
# called.samples <- colSums(!is.na(dp))
# avg.RD <- sum.RD / called.samples



# ### mel:
OR <- seqGetData(array_genofile, "annotation/format/OR")
SR <- seqGetData(array_genofile, "annotation/format/SR")
sum.RD <- colSums(OR + SR, na.rm = TRUE)
called.samples <- colSums(!is.na(OR + SR))
# avg.RD <- colMeans(OR + SR, na.rm = TRUE)
avg.RD <- sum.RD / called.samples

avg.OR <- colMeans(OR, na.rm = TRUE)
avg.SR <- colMeans(SR, na.rm = TRUE)


####################################################################

### per-site heterozygosity
gt <- seqGetData(array_genofile, "genotype")
print(dim(gt))
a1 <- gt[1, , ]
a2 <- gt[2, , ]
# a1 <- matrix(a1, nrow = nvar)
# a2 <- matrix(a2, nrow = nvar)

# called_gt means both alleles are present
called_gt <- !is.na(a1) & !is.na(a2) & a1 >= 0 & a2 >= 0

# genotype categories 
homo_ref <- called_gt & (a1==0) & (a2==0)
homo_alt <- called_gt & (a1==1) & (a2==1)
het_gt <- called_gt & (a1 != a2)

n_called_gt <- colSums(called_gt)
n_homo_ref <- colSums(homo_ref)
n_homo_alt <- colSums(homo_alt)
n_het <- colSums(het_gt)

prop_het <- ifelse(
    n_called_gt > 0,
    n_het / n_called_gt,
    NA_real_
)

### recalculate AF 
af = seqAlleleFreq(array_genofile, ref.allele = 1L)


### Add FIS (ALAN)
# FIS is like FST but for individual genotypes, comparing observed heterozygosity to what is expected under H-W

# p is allele frequency for the snp 
fis <- function(nHet, nTot, p) {
    numerator <- (nHet/nTot)
    denominator <- (nTot/(nTot-1))*(2*p*(1-p) - (nHet/nTot)/(2*nTot))
    1 - numerator/denominator
  }

site_fis <- fis(
    nHet = n_het,
    nTot = n_called_gt,
    p = af
)

chr <- seqGetData(array_genofile, "chromosome")
pos <- seqGetData(array_genofile, "position")

site_rd <- data.frame(
    chr = chr,
    pos = pos, 
    variant.id = variants,
    ### mel only:
    avg.OR = avg.OR,
    avg.SR = avg.SR,
    ### end mel only ^
    avg.RD = avg.RD,
    sum.RD = sum.RD,
    called.samples = called.samples,
    af = af,
    n_homo_ref=n_homo_ref, # ref/ref
    n_het = n_het, # number of samples with heterozygous genotype
    n_homo_alt=n_homo_alt, # alt/alt
    n_called_gt = n_called_gt, # number of samples with non-missing genotype
    prop_het = prop_het, # observed heterozygosity
    fis = site_fis
)

head(site_rd)
summary(site_rd$avg.RD)
summary(site_rd$sum.RD)
summary(site_rd$fis)
summary(site_rd$prop_het)


saveRDS(site_rd, outfile)
seqClose(array_genofile)


# # ####################################
# # # combine intermediate files - do outside of loop
# library(data.table)

# outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"
# files <- list.files(
#     outdir,
#     pattern = "^sim_site_RD_[0-9]+_[0-9]+\\.rds$",
#     full.names = TRUE
# )

# files <- sort(files)
# print(files)
# site_rd_list <- lapply(files, readRDS)

# site_rd <- rbindlist(site_rd_list)

# saveRDS(site_rd, file.path(outdir, "sim_site_RD.rds"))

# print(site_rd)
# print(dim(site_rd))