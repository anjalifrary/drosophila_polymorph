library(SeqArray)

vcf_file <- "/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2.source_BCM-HGSC.dm6.final.norm.ann.eff.vcf.gz"
# gds_file <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/DGRP2.source_BCM-HGSC.dm6.final.reheadered.ann.eff.gds"
gds_file <- "/scratch/ejy4bu/drosophila/gds_files/DGRP2.source_BCM-HGSC.dm6.final.norm.ann.eff.gds"

seqVCF2GDS(vcf.fn = vcf_file, out.fn = gds_file, parallel = 8)

#open and check file
gds <- seqOpen(gds_file)
print(seqSummary(gds))
seqClose(gds)

file.copy(from = gds_file, to = "/project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/")
