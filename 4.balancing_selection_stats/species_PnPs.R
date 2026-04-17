library(SeqArray)
library(data.table)

mel_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/species_rdsFiles/mel_filtered_eff_snp_dt.rds")
sim_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/species_rdsFiles/sim_filtered_eff_snp_dt.rds")

mel_syn <- mel_rds[effect == "synonymous_variant"]
mel_nonsyn <- mel_rds[effect == "missense_variant"]
# message(nrow(mel_rds), " total mel variants")
message(nrow(mel_syn), " synonymous mel variants")
message(nrow(mel_nonsyn), " nonsynonymous mel variants")
message((nrow(mel_syn)+nrow(mel_nonsyn)), " = ", (nrow(mel_rds)), " total mel variants")

sim_syn <- sim_rds[effect == "synonymous_variant"]
sim_nonsyn <- sim_rds[effect == "missense_variant"]
message(nrow(sim_syn), " synonymous sim variants")
message(nrow(sim_nonsyn), " nonsynonymous sim variants")
message((nrow(sim_syn)+nrow(sim_nonsyn)), " = ", (nrow(sim_rds)), " total sim variants")

