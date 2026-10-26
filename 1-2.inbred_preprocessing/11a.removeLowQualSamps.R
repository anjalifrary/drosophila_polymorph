library(SeqArray)
library(data.table)

# ### mel
outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"
# mel_gds <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.gds"
mel_gds <- "/scratch/ejy4bu/drosophila/inbred/gds_files/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.gds"
# out_gds <- paste0(outdir, "DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.goodSamps.gds")
out_gds <- paste0(outdir, "/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.goodSamps.gds")

# meta_file <- as.data.table(read.csv("/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/DGRP2.source_BCM-HGSC.dm6.csv"))
genofile <- seqOpen(mel_gds)
meta_file <- as.data.table(read.csv("/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/DGRP2.source_BCM-HGSC.dm6.10_2026.csv"))

### sim
outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"
sim_gds <- "/scratch/ejy4bu/drosophila/inbred/gds_files/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.dm6.sorted.m2.M2.gds"
out_gds <- paste0(outdir, "dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.dm6.sorted.m2.M2.goodSamps.goodSites.gds")
meta_file <-  as.data.table(read.csv("/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/signor.dsim3.sampleFilt.10_2026.csv"))

genofile <- seqOpen(sim_gds)


sample_ids <- seqGetData(genofile, "sample.id")

bad_samps <- meta_file[Recommendation != "", sample.id]
bad_samps[!bad_samps %in% sample_ids] # should be none

good_samps <- setdiff(sample_ids, bad_samps)

length(sample_ids)
length(bad_samps)
length(good_samps)

seqSetFilter(genofile, sample.id = good_samps)

seqExport(genofile,out_gds)
seqClose(genofile)