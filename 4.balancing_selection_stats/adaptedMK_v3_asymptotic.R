library(SeqArray)
library(data.table)
library(doMC)
registerDoMC(16)

# https://github.com/connor122721/SharedPolymorphismsDaphnia/blob/864eda3b218e1bd2880907a651d791165021261e/Figures/Figure_3_Excess_TSPs/Scripts/Neutral_shared_polymorphism_expectation/1.TSP_MasterScript.R#L10
## for inspo


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
        alpha <- 1 - (Ps * SPns) / (Pns * SPs)
        OR <- 1 - alpha 
        data.table(
            alpha = alpha, 
            OR = OR,
            logOR = log(OR)
        )
    } else {
        data.table(
            alpha = NA_real_,
            OR = NA_real_,
            logOR = NA_real_
        )
    }
}

safe_fread <- function(file) {
    if (!file.exists(file) || file.info(file)$size == 0) {
        return(data.table(chr=character(), pos=integer()))
    }
    fread(file, sep="\t", col.names=c("chr","pos"))
}

MAF_def = "polyAF"
# MAF_def = "globalAF"

# background="mel/sim"
background="mel_only"
maf_inputs <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.49)
# maf_inputs <- c(0.005)
results <- rbindlist(
    foreach(maf = maf_inputs, .packages="data.table") %dopar% {
        maf_label = maf * 100
        # select filtered background
        bg_dir <- "/scratch/ejy4bu/drosophila/GO/gowinda/backgroundFiles/"
        cand_dir <- "/scratch/ejy4bu/drosophila/GO/gowinda/candidateFiles/"

        # bg <- fread(paste0(bg_dir, "MAF", maf_label, "filter_", MAF_def, "/bg_speciesSpecific_", maf_label, "_", MAF_def, ".txt"),
        #     sep="\t", col.names=c("chr", "pos"))
        bg <- safe_fread(paste0(bg_dir, "MAF", maf_label, "filter_", MAF_def, "/bg_speciesSpecific_", maf_label, "_", MAF_def, ".txt"))
        # cand_tsp <- fread(paste0(cand_dir, "MAF", maf_label, "filter_", MAF_def, "/candidate_chrpos_AB_", maf_label, "_", MAF_def, ".txt"),
        #     sep="\t", col.names=c("chr", "pos"))
        cand_tsp <- safe_fread(paste0(cand_dir, "MAF", maf_label, "filter_", MAF_def, "/candidate_chrpos_AB_", maf_label, "_", MAF_def, ".txt"))
        # cand_conv <- fread(paste0(cand_dir, "MAF", maf_label, "filter_", MAF_def, "/candidate_chrpos_FGOPXY_", maf_label, "_", MAF_def, ".txt"),
        #     sep="\t", col.names=c("chr", "pos"))
        cand_conv <- safe_fread(paste0(cand_dir, "MAF", maf_label, "filter_", MAF_def, "/candidate_chrpos_FGOPXY_", maf_label, "_", MAF_def, ".txt"))
        # cand_all <- fread(paste0(cand_dir, "MAF", maf_label, "filter_", MAF_def, "/candidate_chrpos_ABFGOPXY_", maf_label, "_", MAF_def, ".txt"),
        #     sep="\t", col.names=c("chr", "pos")) 
        cand_all <- safe_fread(paste0(cand_dir, "MAF", maf_label, "filter_", MAF_def, "/candidate_chrpos_ABFGOPXY_", maf_label, "_", MAF_def, ".txt"))

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

results <- as.data.table(results)

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
        "OR",
        "logOR"
    )
)


saveRDS(results, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/adaptedMK/new_asymptotic_MK_longResults_polyAF_speciesSpecificBG.rds")

# FIGURES 



# # candidate files of chr + pos (txt) at /scratch/ejy4bu/drosophila/GO/gowinda/candidateFiles/MAFxxfilter_polyAF/candidate_chrpos_AB_xx_polyAF.txt
# # bg files of chr + pos (.txt) at /scratch/ejy4bu/drosophila/GO/gowinda/backgroundFiles/MAFxxfilter_polyAF/bg_speciesSpecific_xx_polyAF.txt

# # use bg and candidate files like the ones made for gowinda (but don't save as text files)
# # add gene info for per-gene stats 

# # run at varying MAF 
# # how to save it? 