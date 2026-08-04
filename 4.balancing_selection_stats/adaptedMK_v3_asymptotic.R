library(SeqArray)
library(data.table)
library(doMC)
registerDoMC(16)

maf_inputs <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.49)

masterCandidates <- readRDS("/scratch/ejy4bu/drosophila/GO/gowinda/candidateFiles/masterCandidateFile.rds")
masterBG_speciesSpecific <- fread("/scratch/ejy4bu/drosophila/GO/gowinda/backgroundFiles/noMAFfilter/bg_speciesSpecific_noMAF.txt",
    sep="\t", col.names=c("chr", "pos"))

master_snp <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/noMAFfilter/all_quality_variants_merge_unfilt.rds")
classed_snp <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/noMAFfilter/subset_qualVar_ofInterest_7-20-2026.rds")
all_snp <- merge(master_snp, classed_snp[, .(chr, pos, classification)], by=c("chr", "pos"), all.x=T)
setcolorder(all_snp, c("chr", "pos", "classification", setdiff(names(all_snp), c("chr", "pos", "classification"))))

masterBG <- merge(masterBG_speciesSpecific, all_snp, by=c("chr", "pos"), all.x=T)

# readRDS("/scratch/ejy4bu/drosophila/GO/gowinda/backgroundFiles/noMAFfilter/bg_speciesSpecific_noMAF.txt")
# masterBG_sharedOnly <- readRDS("/scratch/ejy4bu/drosophila/GO/gowinda/backgroundFiles/noMAFfilter/bg_sharedOnly_noMAF.txt")
# make a BG of mel only + scaled MAF for mel bg 

tsp <- c("A", "B")
conv <- c("F", "G", "O", "P", "X", "Y")

# bg_rds <- masterBG_speciesSpecific

# setDT(bg_rds)
# setindex(bg_rds, gene_id_mel)
# genes <- na.omit(unique(bg_rds$gene_id_mel))
# test subset
# genes <- genes[1:50]

asymptotic_MKlike_stats <- function(Ps, Pns, SPs, SPns, pseudo=0, min_count=0) {
    Ps <- Ps + pseudo
    Pns <- Pns + pseudo
    SPs <- SPs + pseudo
    SPns <- SPns + pseudo

     if (Ps > min_count && Pns > min_count && SPs > min_count && SPns > min_count) {
        OR <- (Pns * SPs) / (Ps * SPns)
        data.table(
            alpha = 1 - 1/OR,
            logOR = log(OR)
        )
    } else {
        data.table(
            alpha = NA_real_,
            logOR = NA_real_
        )
    }
}

MAF_def = "polyAF"
# MAF_def = "globalAF"

# background="mel/sim"
background="mel_only"

results <- rbindlist(
    foreach(maf = maf_inputs, .packages="data.table") %dopar% {
        maf_label = maf * 100
        # select filtered background
        bg_dir <- "/scratch/ejy4bu/drosophila/GO/gowinda/backgroundFiles/"
        cand_dir <- "/scratch/ejy4bu/drosophila/GO/gowinda/candidateFiles/"

        bg <- fread(paste0(bg_dir, "MAF", maf_label, "filter_", MAF_def, "/bg_speciesSpecific_", maf_label, "_", MAF_def, ".txt"),
            sep="\t", col.names=c("chr", "pos"))
        cand_tsp <- fread(paste0(cand_dir, "MAF", maf_label, "filter_", MAF_def, "/candidate_chrpos_AB_", maf_label, "_", MAF_def, ".txt"),
            sep="\t", col.names=c("chr", "pos"))
        cand_conv <- fread(paste0(cand_dir, "MAF", maf_label, "filter_", MAF_def, "/candidate_chrpos_FGOPXY_", maf_label, "_", MAF_def, ".txt"),
            sep="\t", col.names=c("chr", "pos"))
        cand_all <- fread(paste0(cand_dir, "MAF", maf_label, "filter_", MAF_def, "/candidate_chrpos_ABFGOPXY_", maf_label, "_", MAF_def, ".txt"),
            sep="\t", col.names=c("chr", "pos"))            

        bg_table <- merge(bg, masterBG, by=c("chr", "pos"), all.x=TRUE)

        tsp_table <- merge(cand_tsp, masterCandidates, by=c("chr","pos"), all.x=TRUE)
        conv_table <- merge(cand_conv, masterCandidates, by=c("chr","pos"), all.x=TRUE)
        both_table <- merge(cand_all, masterCandidates, by=c("chr","pos"), all.x=TRUE)

        setindex(bg_table, gene_id_mel)
        genes <- unique(na.omit(bg_table$gene_id_mel))

        maf_results <- rbindlist(lapply(genes, function(gene) {
            bg_dt <- bg_table[gene_id_mel == gene]

            Ps <- nrow(bg_dt[!is.na(ref_mel) & is.na(ref_sim) & effect_mel %like% "syn" ])
            Pns <- nrow(bg_dt[!is.na(ref_mel) & is.na(ref_sim) & effect_mel %like% "missense"])

            get_SP <- function(candidate_dt){
                cand_dt <- candidate_dt[gene_id_mel == gene]

                SPs <- nrow(cand_dt[effect_mel %like% "syn"])
                SPns <- nrow(cand_dt[effect_mel %like% "missense"])

                asymptotic_MKlike_stats(Ps, Pns, SPs, SPns)[, ':=' (
                    Ps = Ps,
                    Pns = Pns, 
                    SPs = SPs, 
                    SPns = SPns
                )] 
            }
            rbind(
            get_SP(tsp_table)[, ':=' (
                    gene=gene, 
                    MAF=maf, 
                    group="TSP",
                    background=background,
                    MAF_def=MAF_def
                )],
            get_SP(conv_table)[, ':=' (
                    gene=gene, 
                    MAF=maf, 
                    group="CONV",
                    background=background,
                    MAF_def=MAF_def
                )],
            get_SP(both_table)[, ':=' (
                    gene=gene, 
                    MAF=maf, 
                    group="BOTH",
                    background=background,
                    MAF_def=MAF_def
                )]
            )

        })
        )
        maf_results

    }
)

setcolorder(
    results,
    c(
        "gene",
        "MAF",
        "group",
        "background",
        "MAF_def",
        "Ps",
        "Pns",
        "SPs",
        "SPns",
        "alpha",
        "logOR"
    )
)


saveRDS(results, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/adaptedMK/asymptotic_MK_longResults_polyAF_speciesSpecificBG.rds")

# FIGURES 





# alpha_beta <- function(Ps, Pns, SPs, SPns, min_count = 0) {
#     pseudo <- 1

#     # Ps <- Ps + pseudo
#     # Pns <- Pns + pseudo
#     SPs <- SPs + pseudo
#     SPns <- SPns + pseudo
#     if (Pns > min_count && SPs > min_count && SPns > min_count && Ps > min_count) 
#     return(1 - (Ps * SPns) / (Pns * SPs))
#     NA_real_
# }


# count_class <- function(sub_dt, effect_pattern, xy_classes = c("X", "Y")) {
#     n_normal <- nrow(sub_dt[!classification %in% xy_classes & effect_sim %like% effect_pattern]) 
#     n_xy     <- nrow(sub_dt[ classification %in% xy_classes & effect_sim %like% effect_pattern]) # /2 stopped dividing by 2 because only 1 of the two rows has a effect_mel; the other has effect_sim
#     n_normal + n_xy
# }

# results <- rbindlist(
#     foreach(gene = genes, .packages = "data.table") %dopar% {

#         dt <- bg_rds[gene_id_mel == gene]
        
#         # # mel species-specific background (same for all tests)
#         mel_Ps  <- dt[!is.na(ref_mel) & is.na(ref_sim) & effect_mel %like% "syn"]
#         mel_Pns <- dt[!is.na(ref_mel) & is.na(ref_sim) & effect_mel %like% "missense"]

#         # sim species-specific background 

#         # sim-only rows don't have gene_id_mel so added logic to define mel genespace for sim Pn/Ps
#         # sim_dt <- bg_rds[
#         #     !is.na(ref_sim) & is.na(ref_mel) &
#         #     chr == unique(dt$chr) &
#         #     pos >= min(dt$pos) &
#         #     pos <= max(dt$pos)
#         # ]

#         # sim_Ps <- sim_dt[effect_sim %like% "syn"]
#         # sim_Pns <- sim_dt[effect_sim %like% "missense"]

#         # Ps  <- nrow(sim_Ps) 
#         # Pns <- nrow(sim_Pns) 

#         # TSP candidates (A+B)
#         tsp_SPs   <- count_class(dt[classification %in% tsp],        "syn")
#         tsp_SPns  <- count_class(dt[classification %in% tsp],        "missense")
#         alpha_tsp <- alpha_beta(Ps, Pns, tsp_SPs, tsp_SPns)

#         # convergent (F,G,O,P,X,Y) 
#         conv_SPs  <- count_class(dt[classification %in% conv],       "syn")
#         conv_SPns <- count_class(dt[classification %in% conv],       "missense")
#         alpha_conv <- alpha_beta(Ps, Pns, conv_SPs, conv_SPns)

#         # all classified together (A+B+F+G+O+P+X+Y)
#         all_SPs   <- count_class(dt[classification %in% c(tsp,conv)],"syn")
#         all_SPns  <- count_class(dt[classification %in% c(tsp,conv)],"missense")
#         alpha_all <- alpha_beta(Ps, Pns, all_SPs, all_SPns)

#         data.table(
#             gene       = gene,
#             alpha_tsp  = alpha_tsp,
#             alpha_conv = alpha_conv,
#             alpha_all  = alpha_all,
#             mel_Ps     = Ps,
#             mel_Pns    = Pns,
#             # sim_Ps     = Ps,
#             # sim_Pns    = Pns,
#             tsp_SPs    = tsp_SPs,
#             tsp_SPns   = tsp_SPns,
#             conv_SPs   = conv_SPs,
#             conv_SPns  = conv_SPns,
#             all_SPs    = all_SPs,
#             all_SPns   = all_SPns
#         )
#     }
# )

# results_clean <- results[!is.na(alpha_tsp) | !is.na(alpha_conv) | !is.na(alpha_all)]

# summary(results_clean$alpha_tsp)
# summary(results_clean$alpha_conv)
# summary(results_clean$alpha_all)

# table(results$tsp_SPns)
# table(results$tsp_SPs)
# table(results$conv_SPns)
# table(results$conv_SPs)


# # saveRDS(results_clean, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results_sim_noPseudoCount_06-29-2026.rds")
# # saveRDS(results_clean, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results_sim_1PseudoAll_06-29-2026.rds")
# # saveRDS(results_clean, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results_sim_1PseudoSPxOnly_06-29-2026.rds")

# ##############################################################################
# ### Get Gene List ### 
# ##############################################################################

# gene_list <- all_df[alpha>0, .(gene, category, alpha, background)]
# gene_list <- gene_list[category != "ALL"]
# gene_wide <- dcast(
#     gene_list,
#     gene ~ background + category,
#     value.var = "alpha"
# )


# # candidate files of chr + pos (txt) at /scratch/ejy4bu/drosophila/GO/gowinda/candidateFiles/MAFxxfilter_polyAF/candidate_chrpos_AB_xx_polyAF.txt
# # bg files of chr + pos (.txt) at /scratch/ejy4bu/drosophila/GO/gowinda/backgroundFiles/MAFxxfilter_polyAF/bg_speciesSpecific_xx_polyAF.txt

# # use bg and candidate files like the ones made for gowinda (but don't save as text files)
# # add gene info for per-gene stats 

# # run at varying MAF 
# # how to save it? 