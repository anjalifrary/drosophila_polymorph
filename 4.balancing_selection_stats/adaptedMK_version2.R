library(SeqArray)
library(data.table)
library(doMC)
registerDoMC(16)


var_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/merged_tables/quality/all_quality_variants_merge_unfilt.rds")
nrow(var_rds)

setDT(var_rds)
setindex(var_rds, gene_id_mel)
genes <- na.omit(unique(var_rds$gene_id_mel))
# test subset
genes <- genes[1:50]

# saveRDS(genes, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_listOfGenes.rds")

# genes <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_listOfGenes.rds")

# gene <- genes[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
# dt <- var_rds[gene_id_mel == gene]

results <- rbindlist(
    foreach(gene=genes, .packages="data.table") 
    %dopar%{
# rbindlist(lapply(genes, function(gene){

    dt <- var_rds[.(gene), on = "gene_id_mel"]

    mel_Ps <- dt[!is.na(ref_mel) & is.na(ref_sim) & effect_mel%like%c("syn")]
    mel_Pns <- dt[!is.na(ref_mel) & is.na(ref_sim) & effect_mel%like%c("missense")]

    sim_Ps <- dt[!is.na(ref_sim) & is.na(ref_mel) & effect_sim%like%c("syn")]
    sim_Pns <- dt[!is.na(ref_sim) & is.na(ref_mel) & effect_sim%like%c("missense")]
    # message("mel Ps = ", nrow(mel_Ps))
    # message("mel Pns = ", nrow(mel_Pns))
    # message("sim Ps = ", nrow(sim_Ps))
    # message("sim Pns = ", nrow(sim_Pns))

    shared_SPs <- dt[!is.na(ref_mel) & !is.na(ref_sim) & effect_mel%like%c("syn") & effect_sim%like%c("syn")]
    shared_SPns <- dt[!is.na(ref_mel) & !is.na(ref_sim) & effect_mel%like%c("missense") & effect_sim%like%c("missense")]
    # message("SPs = ", nrow(shared_SPs))
    # message("SPns = ", nrow(shared_SPns))

    SPs <- nrow(shared_SPs)
    SPns <- nrow(shared_SPns)

    # mel alpha beta 
    Ps <- nrow(mel_Ps)
    Pns <- nrow(mel_Pns)
    alpha_beta_mel <- NA_real_
    if (Pns > 0 & SPs > 0) {alpha_beta_mel <- 1 - (Ps * SPns)/(Pns * SPs)}

    # sim alpha beta
    Ps <- nrow(sim_Ps)
    Pns <- nrow(sim_Pns)
    alpha_beta_sim <- NA_real_
    if (Pns > 0 & SPs > 0) {alpha_beta_sim <- 1 - (Ps * SPns)/(Pns * SPs)}

    data.table(
        gene = gene,
        mel_Ps = nrow(mel_Ps),
        mel_Pns = nrow(mel_Pns),
        sim_Ps = nrow(sim_Ps),
        sim_Pns = nrow(sim_Pns),
        SPs = SPs,
        SPns = SPns,
        alpha_beta_mel = alpha_beta_mel,
        alpha_beta_sim = alpha_beta_sim
    )
}
)
# ))

saveRDS(
    results,
    "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results.rds"
)

fwrite(
    results,
    "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results.txt",
    sep="\t"
)