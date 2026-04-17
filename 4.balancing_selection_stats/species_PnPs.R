library(SeqArray)
library(data.table)

mel_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/species_rdsFiles/mel_filtered_eff_snp_dt.rds")
sim_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/species_rdsFiles/sim_filtered_eff_snp_dt.rds")
merged_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variantsOfInterest_qualFiltered.rds")
# NOTE: I am using the variants of subset datatable which is already quality filtered for both species, includes all variants at same site and variants within codon
# but check if same results are obtained when calculating from recreated merged datatable 

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

# setDT(merged_rds)

shared_variants <- merged_rds[
    !is.na(ref_mel) & 
    !is.na(ref_sim)
]

shared_syn <- shared_variants[
    (alt_mel == alt_sim)
]

shared_nonsyn <- shared_variants[
    (alt_mel != alt_sim)
]

message(nrow(merged_rds), " variants in merged datatable")
message(nrow(shared_variants), " shared variants")
message(nrow(shared_syn), " synonymous shared variants")
message(nrow(shared_nonsyn), " nonsynonymous shared variants")
message((nrow(shared_syn)+nrow(shared_nonsyn)), " = ", (nrow(shared_variants)), " total shared variants")