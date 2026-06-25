library(SeqArray)
library(data.table)
library(doMC)
registerDoMC(16)


var_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/merged_tables/quality/all_quality_variants_merge_unfilt.rds")
nrow(var_rds)

classed_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed.rds")
nrow(classed_rds)
names(classed_rds)

tsp <- c("A", "B")
conv <- c("F", "G", "O", "P", "X", "Y")

var_rds <- merge(var_rds, classed_rds[, .(chr, pos, classification)], by=c("chr", "pos"), all.x=T)
nrow(var_rds[is.na(classification)])
nrow(var_rds[!is.na(classification)])


setDT(var_rds)
setindex(var_rds, gene_id_mel)
genes <- na.omit(unique(var_rds$gene_id_mel))
# test subset
# genes <- genes[1:50]


## what should the min count be? this is crucial for the convergent test bc lots of 0s for SPns
alpha_beta <- function(Ps, Pns, SPs, SPns, min_count = 0) {
    if (Pns > min_count && SPs > min_count) return(1 - (Ps * SPns) / (Pns * SPs))
    NA_real_
}

count_class <- function(sub_dt, effect_pattern, xy_classes = c("X", "Y")) {
    n_normal <- nrow(sub_dt[!classification %in% xy_classes & effect_mel %like% effect_pattern])
    n_xy     <- nrow(sub_dt[ classification %in% xy_classes & effect_mel %like% effect_pattern]) / 2
    n_normal + n_xy
}

results <- rbindlist(
    foreach(gene = genes, .packages = "data.table") %dopar% {

        dt <- var_rds[.(gene), on = "gene_id_mel"]

        # mel species-specific background (same for all tests)
        mel_Ps  <- dt[!is.na(ref_mel) & is.na(ref_sim) & effect_mel %like% "syn"]
        mel_Pns <- dt[!is.na(ref_mel) & is.na(ref_sim) & effect_mel %like% "missense"]
        Ps  <- nrow(mel_Ps)
        Pns <- nrow(mel_Pns)

        # TSP candidates (A+B)
        tsp_SPs   <- count_class(dt[classification %in% tsp],        "syn")
        tsp_SPns  <- count_class(dt[classification %in% tsp],        "missense")
        alpha_tsp <- alpha_beta(Ps, Pns, tsp_SPs, tsp_SPns)

        # convergent (F,G,O,P,X,Y) 
        conv_SPs  <- count_class(dt[classification %in% conv],       "syn")
        conv_SPns <- count_class(dt[classification %in% conv],       "missense")
        alpha_conv <- alpha_beta(Ps, Pns, conv_SPs, conv_SPns)

        # all classified together (A+B+F+G+O+P+X+Y)
        all_SPs   <- count_class(dt[classification %in% c(tsp,conv)],"syn")
        all_SPns  <- count_class(dt[classification %in% c(tsp,conv)],"missense")
        alpha_all <- alpha_beta(Ps, Pns, all_SPs, all_SPns)

        data.table(
            gene       = gene,
            alpha_tsp  = alpha_tsp,
            alpha_conv = alpha_conv,
            alpha_all  = alpha_all,
            mel_Ps     = Ps,
            mel_Pns    = Pns,
            tsp_SPs    = tsp_SPs,
            tsp_SPns   = tsp_SPns,
            conv_SPs   = conv_SPs,
            conv_SPns  = conv_SPns,
            all_SPs    = all_SPs,
            all_SPns   = all_SPns
        )
    }
)

results <- results[!is.na(alpha_tsp) | !is.na(alpha_conv) | !is.na(alpha_all)]

saveRDS(results, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results.rds")
# fwrite(results, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results.txt", sep="\t")
