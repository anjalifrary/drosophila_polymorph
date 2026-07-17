library(SeqArray)

vcf_file <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/DGRP2.source_BCM-HGSC.dm6.final.ann.eff.vcf.gz"
gds_file <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/DGRP2.source_BCM-HGSC.dm6.final.ann.eff.gds"

seqVCF2GDS(vcf.fn = vcf_file, out.fn = gds_file, parallel = TRUE)

#open and check file
gds <- seqOpen(gds_file)
print(seqSummary(gds))
seqClose(gds)