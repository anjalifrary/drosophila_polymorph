library(SeqArray)

# ## sim:
# vcf_file <- "/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.vcf.gz"
# gds_file <- "/scratch/ejy4bu/drosophila/inbred/gds_files/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.gds"

## mel:
# vcf_file <- "/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.vcf.gz"
# gds_file <- "/scratch/ejy4bu/drosophila/inbred/gds_files/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.gds"

seqVCF2GDS(vcf.fn = vcf_file, out.fn = gds_file, parallel = 8)

#open and check file
gds <- seqOpen(gds_file)
print(seqSummary(gds))
seqClose(gds)

file.copy(from = gds_file, to = "/project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/")
