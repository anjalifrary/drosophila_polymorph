library(SeqArray)

vcf_file <- "/scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.fixed.vcf.gz"
gds_file <- "/scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds"

seqVCF2GDS(vcf.fn = vcf_file, out.fn = gds_file, parallel = TRUE)

#open and check file
gds <- seqOpen(gds_file)
print(seqSummary(gds))
seqClose(gds)