library(SeqArray)
library(data.table)

#################################################################################################
# This script takes a datatable (rds) of variants, applies a quality filter, and recalculates AF 
#################################################################################################

# load original gds files
mel_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.eff.gds")
sim_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.eff.gds")

# load unfiltered merged rds 
unfilt_dt <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/raw_merged_tables/all_variants_merge_unfilt.rds")

# load metadata files
mel_metadata <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/dest_v2.samps_24Aug2024.xa.csv")
# sim_metadata <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/simulans_pooled.meta.use.csv")
sim_contam <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/Dsim_contamRates.csv")

# load variants of interest rds file
filtered_dt <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variantsOfInterest.rds")


######## Get variant ids for subset of variants chosen ##################

setkey(unfilt_dt, chr, pos) # join on chr + pos 
setkey(filtered_dt, chr, pos)

variants_of_interest <- unfilt_dt[filtered_dt, on = .(chr, pos), nomatch = 0L]

mel_var_ids <- variants_of_interest[!is.na(variant.id_mel), variant.id_mel]
sim_var_ids <- variants_of_interest[!is.na(variant.id_sim), variant.id_sim]

# mel_var_ids <- na.omit(variants_of_interest$variant.id_mel)
# sim_var_ids <- na.omit(variants_of_interest$variant.id_sim)

message(length(mel_var_ids), " mel variant IDs to extract")
message(length(sim_var_ids), " sim variant IDs to extract")


#################################################################################################
# 1. Filter melanogaster samples, by metadata recommendation 

### filter mel samps
good_mel_samps <- mel_metadata[Recommendation == "Pass", sampleId]
good_mel_samps <- intersect(seqGetData(mel_gds, "sample.id"), good_mel_samps)
message(length(good_mel_samps), " mel samples pass filter")

seqSetFilter(mel_gds, sample.id = good_mel_samps, variant.id = mel_var_ids)

#################################################################################################
# 2. Filter simulans samples, by contamination rate (5% contamination limit) 

good_sim_samps <- sim_contam[propSim >= 0.95, sampleId]
good_sim_samps <- intersect(seqGetData(sim_gds, "sample.id"), good_sim_samps)
seqSetFilter(sim_gds, sample.id = good_sim_samps, variant.id = sim_var_ids)

message(length(good_sim_samps), " sim samples pass filter")

#################################################################################################
# 3. Recalculate allele frequency

mel_af_dt <- data.table(
  variant.id_mel = seqGetData(mel_gds, "variant.id"),
  af_mel_recomputed = 1 - seqAlleleFreq(mel_gds) # had to flip because it was calculating the allele frequency of the reference? 
)

sim_af_dt <- data.table(
  variant.id_sim = seqGetData(sim_gds, "variant.id"),
  af_sim_recomputed = seqAlleleFreq(sim_gds)
)

##############################################################################
# 4. Rejoin tables by qual filter samples

# join recalculated AF back to variants of interest by variant ID
variants_of_interest <- mel_af_dt[variants_of_interest, on = .(variant.id_mel)]
variants_of_interest <- sim_af_dt[variants_of_interest, on = .(variant.id_sim)]

#### Check allele frequencies
variants_of_interest[!is.na(af_mel) & (af_mel_recomputed>0), .(
  mean_af_mel_before = mean(af_mel), mean_af_mel_after = mean(af_mel_recomputed),
  mean_diff = mean(abs(af_mel - af_mel_recomputed))
)]

variants_of_interest[!is.na(af_sim) & (af_sim_recomputed>0), .(
  mean_af_sim_before = mean(af_sim), mean_af_sim_after = mean(af_sim_recomputed),
  mean_diff = mean(abs(af_sim - af_sim_recomputed))
)]

# filter! 
quality_dt <- variants_of_interest[
  (is.na(ref_mel) | (af_mel_recomputed>0)) &
  (is.na(ref_sim) | (af_sim_recomputed>0))
]
message(nrow(variants_of_interest), " variants before quality filter")
message(nrow(quality_dt), " variants after quality filter")

#############################################################################
# 5. Save file


quality_dt[, af_mel := af_mel_recomputed]
quality_dt[, af_sim := af_sim_recomputed]
quality_dt[, af_mel_recomputed := NULL]
quality_dt[, af_sim_recomputed := NULL]

########## Fix the column names 
cols <- names(filtered_dt)
# setcolorder(quality_dt, cols)

quality_dt <- quality_dt[, ..cols]
names(quality_dt)

# cols_to_remove <- grep("^i\\.", names(quality_dt), value = TRUE)
# quality_dt[, (cols_to_remove) := NULL]
# new_cols <- names(quality_dt)


out_csv <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variantsOfInterest_qualFiltered_test500.csv"
subset_table <- quality_dt[1:500, ]
fwrite(subset_table, out_csv)
message("saved first 500 rows to csv at ", out_csv)

saveRDS(quality_dt, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variantsOfInterest_qualFiltered.rds")

seqClose(mel_gds)
seqClose(sim_gds)
message("complete!")
