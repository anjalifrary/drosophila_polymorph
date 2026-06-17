library(SeqArray)
library(data.table)

# see adaptedMK_version2.R for more updated version

# load species-specific tables
# NOTE: I am using the variants of subset datatable which is already quality filtered for both species, includes all variants at same site and variants within codon
    # but check if same results are obtained when calculating from recreated merged datatable   
mel_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/species_rdsFiles/mel_filtered_eff_snp_dt.rds")
sim_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/species_rdsFiles/sim_filtered_eff_snp_dt.rds")
merged_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variantsOfInterest_qualFiltered.rds")

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




############################################################################3
# redoing: 6/17/2026
var_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/merged_tables/quality/all_quality_variants_merge_unfilt.rds")
nrow(var_rds)

# genes <- na.omit(unique(var_rds$gene_id_mel))
# saveRDS(genes, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_listOfGenes.rds")

genes <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_listOfGenes.rds")

gene <- genes[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
dt <- var_rds[gene_id_mel == gene]

mel_Ps <- dt[!is.na(ref_mel) & is.na(ref_sim) & effect_mel%like%c("syn")]
mel_Pns <- dt[!is.na(ref_mel) & is.na(ref_sim) & effect_mel%like%c("missense")]

sim_Ps <- dt[!is.na(ref_sim) & is.na(ref_mel) & effect_sim%like%c("syn")]
sim_Pns <- dt[!is.na(ref_sim) & is.na(ref_mel) & effect_sim%like%c("missense")]
message("mel Ps = ", nrow(mel_Ps))
message("mel Pns = ", nrow(mel_Pns))
message("sim Ps = ", nrow(sim_Ps))
message("sim Pns = ", nrow(sim_Pns))

shared_SPs <- dt[!is.na(ref_mel) & !is.na(ref_sim) & effect_mel%like%c("syn") & effect_sim%like%c("syn")]
shared_SPns <- dt[!is.na(ref_mel) & !is.na(ref_sim) & effect_mel%like%c("missense") & effect_sim%like%c("missense")]
message("SPs = ", nrow(shared_SPs))
message("SPns = ", nrow(shared_SPns))

Ps <- nrow(mel_Ps)
Pns <- nrow(mel_Pns)
SPs <- nrow(shared_SPs)
SPns <- nrow(shared_SPns)

alpha_beta_mel <- NA_real_
if (Pns > 0 & SPs > 0) {alpha_beta_mel <- 1 - (Ps * SPns)/(Pns * SPs)}

Ps <- nrow(sim_Ps)
Pns <- nrow(sim_Pns)
alpha_beta_sim <- NA_real_
if (Pns > 0 & SPs > 0) {alpha_beta_sim <- 1 - (Ps * SPns)/(Pns * SPs)}

# moving this stuff to adaptedMK_version2
############################# 


# genome-wide summary
mel_Ps <- var_rds[!is.na(ref_mel) & is.na(ref_sim) & effect_mel%like%c("syn")]
mel_Pns <- var_rds[!is.na(ref_mel) & is.na(ref_sim) & effect_mel%like%c("missense")]

sim_Ps <- var_rds[!is.na(ref_sim) & is.na(ref_mel) & effect_sim%like%c("syn")]
sim_Pns <- var_rds[!is.na(ref_sim) & is.na(ref_mel) & effect_sim%like%c("missense")]
message("mel Ps = ", nrow(mel_Ps))
message("mel Pns = ", nrow(mel_Pns))
message("sim Ps = ", nrow(sim_Ps))
message("sim Pns = ", nrow(sim_Pns))

shared_SPs <- var_rds[!is.na(ref_mel) & !is.na(ref_sim) & effect_mel%like%c("syn") & effect_sim%like%c("syn")]
shared_SPns <- var_rds[!is.na(ref_mel) & !is.na(ref_sim) & effect_mel%like%c("missense") & effect_sim%like%c("missense")]
message("SPs = ", nrow(shared_SPs))
message("SPns = ", nrow(shared_SPns))

# which group do these cases fall under? 
nrow(var_rds[!is.na(ref_mel) & !is.na(ref_sim) & effect_mel%like%c("missense") & effect_sim%like%c("syn")])
nrow(var_rds[!is.na(ref_mel) & !is.na(ref_sim) & effect_mel%like%c("syn") & effect_sim%like%c("missense")])

nrow(var_rds)
sum(nrow(mel_Ps), nrow(mel_Pns), nrow(sim_Ps), nrow(sim_Pns), nrow(shared_SPs), nrow(shared_SPns))

num <- prod(nrow(mel_Ps), nrow(shared_SPns))
denom <- prod(nrow(mel_Pns), nrow(shared_SPs))
alpha_beta <- 1 - (num)/(denom)
message(alpha_beta)
# if alpha_beta > 0, indicates balancing selection