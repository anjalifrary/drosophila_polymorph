library(SeqArray)
library(data.table)

# load gds files
mel_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.eff.gds")
sim_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.eff.gds")

# load metadata files
mel_metadata <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/dest_v2.samps_24Aug2024.xa.csv")
sim_metadata <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/simulans_pooled.meta.use.csv")
sim_contam <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/Dsim_contamRates.csv")


### filter mel samps
bad_mel_samps <- mel_metadata[Recommendation != "Pass", sampleId]

# seqGetData(mel_gds, "sample.id")[1:10]

# get the good samples
good_mel_samps <- setdiff(seqGetData(mel_gds, "sample.id"), bad_mel_samps)

seqSetFilter(mel_gds, sample.id = good_mel_samps)

# recalculate af based on samples that pass filter
af <- seqAlleleFreq(mel_gds)

variant.id <- seqGetData(mel_gds, "variant.id")
chr <- seqGetData(mel_gds, "chromosome")
pos <- seqGetData(mel_gds, "position")

dt <- data.table(
  chr = chr,
  pos = pos,
  variant.id = variant.id,
  af_mel = af
)

### idk what next
