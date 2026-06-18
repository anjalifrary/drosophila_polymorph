library(SeqArray)
library(data.table)

#################################################################################################
# This script takes two rds files for species, applies a quality filter, and recalculates AF for updated species tables
#################################################################################################

# load original gds files
mel_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.eff.gds")
sim_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.eff.gds")

# load rds files NOTE: this is only synonymous and missense. should i use all?
mel_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/species_rdsFiles/mel_eff_snp_dt.rds")
sim_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/species_rdsFiles/sim_eff_snp_dt.rds")

# load metadata files
mel_metadata <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/dest_v2.samps_24Aug2024.xa.csv")
# sim_metadata <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/simulans_pooled.meta.use.csv")
sim_contam <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/Dsim_contamRates.csv")

######## Get variant ids for subset of variants chosen ##################

mel_var_ids <- mel_rds[, variant.id]
sim_var_ids <- sim_rds[, variant.id]

message(length(mel_var_ids), " mel variant IDs to extract")
message(length(sim_var_ids), " sim variant IDs to extract")

#################################################################################################
# 1. Filter melanogaster samples, by metadata recommendation 

### filter mel samps
good_mel_samps <- mel_metadata[Recommendation == "Pass", sampleId]
good_mel_samps <- intersect(seqGetData(mel_gds, "sample.id"), good_mel_samps)
message(length(good_mel_samps), " mel samples pass filter")

length(intersect(mel_var_ids, seqGetData(mel_gds, "variant.id")))
seqSetFilter(mel_gds, sample.id = good_mel_samps, variant.id = mel_var_ids)

#################################################################################################
# 2. Filter simulans samples, by contamination rate (5% contamination limit) 

good_sim_samps <- sim_contam[propSim >= 0.95, sampleId]
good_sim_samps <- intersect(seqGetData(sim_gds, "sample.id"), good_sim_samps)
seqSetFilter(sim_gds, sample.id = good_sim_samps, variant.id = sim_var_ids)

length(intersect(sim_var_ids, seqGetData(sim_gds, "variant.id")))
message(length(good_sim_samps), " sim samples pass filter")

#################################################################################################
# 3. Recalculate allele frequency

mel_af_dt <- data.table(
  variant.id = seqGetData(mel_gds, "variant.id"),
#   af_recomputed = seqAlleleFreq(mel_gds)
  af_recomputed = 1 - seqAlleleFreq(mel_gds) # had to flip because it was calculating the allele frequency of the reference? 
)

sim_af_dt <- data.table(
  variant.id = seqGetData(sim_gds, "variant.id"),
  af_recomputed = seqAlleleFreq(sim_gds)
)

##############################################################################
# 4. Rejoin tables by qual filter samples

# join recalculated AF back to variants of interest by variant ID
mel_rds <- mel_af_dt[mel_rds, on = .(variant.id)]
sim_rds <- sim_af_dt[sim_rds, on = .(variant.id)]

#### Check allele frequencies
mel_rds[!is.na(af) & (af_recomputed>0), .(
  mean_af_mel_before = mean(af), mean_af_mel_after = mean(af_recomputed),
  mean_diff = mean(abs(af - af_recomputed))
)]

sim_rds[!is.na(af) & (af_recomputed>0), .(
  mean_af_sim_before = mean(af), mean_af_sim_after = mean(af_recomputed),
  mean_diff = mean(abs(af - af_recomputed))
)]

# filter! 
mel_new_dt <- mel_rds[
    (af_recomputed>0) & (af_recomputed<1)
]
message(nrow(mel_rds), " mel variants BEFORE filter")
message(nrow(mel_new_dt), " mel variants AFTER filter")

sim_new_dt <- sim_rds[
    (af_recomputed>0) & (af_recomputed<1)
]
message(nrow(sim_rds), " sim variants BEFORE filter")
message(nrow(sim_new_dt), " sim variants AFTER filter")

#############################################################################
# 5. Save files

mel_new_dt[, af := af_recomputed]
mel_new_dt[, af_recomputed := NULL]

sim_new_dt[, af := af_recomputed]
sim_new_dt[, af_recomputed := NULL]

cols_to_remove <- grep("^i\\.", names(mel_new_dt), value = TRUE)
mel_new_dt[, (cols_to_remove) := NULL]

cols_to_remove <- grep("^i\\.", names(sim_new_dt), value = TRUE)
sim_new_dt[, (cols_to_remove) := NULL]

mel_csv <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/species_rdsFiles/mel_filtered_eff_snp_dt.csv"
subset_mel_table <- mel_new_dt[1:500, ]
fwrite(subset_mel_table, mel_csv)
message("saved first 500 rows to csv at ", mel_csv)

saveRDS(mel_new_dt, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/species_rdsFiles/mel_filtered_eff_snp_dt.rds")


sim_csv <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/species_rdsFiles/sim_filtered_eff_snp_dt.csv"
subset_sim_table <- sim_new_dt[1:500, ]
fwrite(subset_sim_table, sim_csv)
message("saved first 500 rows to csv at ", sim_csv)

saveRDS(sim_new_dt, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/species_rdsFiles/sim_filtered_eff_snp_dt.rds")

seqClose(mel_gds)
seqClose(sim_gds)
message("complete!")
