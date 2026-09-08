library(SeqArray)
library(data.table)

# ### mel
# outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"
# mel_gds <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.gds"
# out_gds <- paste0(outdir, "DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.goodSamps.gds")
# meta_file <- as.data.table(read.csv("/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/DGRP2.source_BCM-HGSC.dm6.csv"))
# genofile <- seqOpen(mel_gds)

### sim
outdir <- "/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/"
sim_gds <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.gds"
out_gds <- paste0(outdir, "dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.goodSamps.gds")
meta_file <-  as.data.table(read.csv("/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/signor.dsim3.sampleFilt.csv"))
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