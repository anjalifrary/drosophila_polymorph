library(SeqArray)
library(data.table)
library(doMC)
registerDoMC(16)


bg_rds <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/all_quality_variants_MAF5_merge_unfilt.rds")
nrow(bg_rds)

classed_rds <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/voi_fromBG_qualVar_ofInterest_MAF5_classed_06-29-2026.rds") # only classes ABFGOPXY
# classed_rds <- readRDS("/project/berglandlab/anjali/drosophila/polymorphism/classification/subset_fromBG_qualVar_ofInterest_MAF5_classed_06-29-2026.rds") # all classes 
nrow(classed_rds)
names(classed_rds)

tsp <- c("A", "B")
conv <- c("F", "G", "O", "P", "X", "Y")
# species <- "sim"

bg_rds <- merge(bg_rds, classed_rds[, .(chr, pos, classification)], by=c("chr", "pos"), all.x=T)
nrow(bg_rds[is.na(classification)])
nrow(bg_rds[!is.na(classification)])


setDT(bg_rds)
setindex(bg_rds, gene_id_mel)
genes <- na.omit(unique(bg_rds$gene_id_mel))
# test subset
# genes <- genes[1:50]


## what should the min count be? this is crucial for the convergent test bc lots of 0s for SPns
alpha_beta <- function(Ps, Pns, SPs, SPns, min_count = 0) {
    pseudo <- 1

    # Ps <- Ps + pseudo
    # Pns <- Pns + pseudo
    SPs <- SPs + pseudo
    SPns <- SPns + pseudo
    if (Pns > min_count && SPs > min_count && SPns > min_count && Ps > min_count) 
    return(1 - (Ps * SPns) / (Pns * SPs))
    NA_real_
}

count_class <- function(sub_dt, effect_pattern, xy_classes = c("X", "Y")) {
    n_normal <- nrow(sub_dt[!classification %in% xy_classes & effect_sim %like% effect_pattern]) 
    n_xy     <- nrow(sub_dt[ classification %in% xy_classes & effect_sim %like% effect_pattern]) # /2 stopped dividing by 2 because only 1 of the two rows has a effect_mel; the other has effect_sim
    n_normal + n_xy
}

results <- rbindlist(
    foreach(gene = genes, .packages = "data.table") %dopar% {

        dt <- bg_rds[gene_id_mel == gene]
        
        # # mel species-specific background (same for all tests)
        # mel_Ps  <- dt[!is.na(ref_mel) & is.na(ref_sim) & effect_mel %like% "syn"]
        # mel_Pns <- dt[!is.na(ref_mel) & is.na(ref_sim) & effect_mel %like% "missense"]

        # sim species-specific background 

        # sim-only rows don't have gene_id_mel so added logic to define mel genespace for sim Pn/Ps
        sim_dt <- bg_rds[
            !is.na(ref_sim) & is.na(ref_mel) &
            chr == unique(dt$chr) &
            pos >= min(dt$pos) &
            pos <= max(dt$pos)
        ]

        sim_Ps <- sim_dt[effect_sim %like% "syn"]
        sim_Pns <- sim_dt[effect_sim %like% "missense"]

        Ps  <- nrow(sim_Ps) 
        Pns <- nrow(sim_Pns) 

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
            # mel_Ps     = Ps,
            # mel_Pns    = Pns,
            sim_Ps     = Ps,
            sim_Pns    = Pns,
            tsp_SPs    = tsp_SPs,
            tsp_SPns   = tsp_SPns,
            conv_SPs   = conv_SPs,
            conv_SPns  = conv_SPns,
            all_SPs    = all_SPs,
            all_SPns   = all_SPns
        )
    }
)

results_clean <- results[!is.na(alpha_tsp) | !is.na(alpha_conv) | !is.na(alpha_all)]

summary(results_clean$alpha_tsp)
summary(results_clean$alpha_conv)
summary(results_clean$alpha_all)

table(results$tsp_SPns)
table(results$tsp_SPs)
table(results$conv_SPns)
table(results$conv_SPs)


# saveRDS(results_clean, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results_sim_noPseudoCount_06-29-2026.rds")
# saveRDS(results_clean, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results_sim_1PseudoAll_06-29-2026.rds")
# saveRDS(results_clean, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results_sim_1PseudoSPxOnly_06-29-2026.rds")

##############################################################################
### Get Gene List ### 
##############################################################################

gene_list <- all_df[alpha>0, .(gene, category, alpha, background)]
gene_list <- gene_list[category != "ALL"]
gene_wide <- dcast(
    gene_list,
    gene ~ background + category,
    value.var = "alpha"
)


##############################################################################
### FIGURES ### 
##############################################################################
library(data.table)
library(ggplot2)

# # no psuedo Bayesian +1 thing
# df_mel <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results_mel_noPseudoCount_06-29-2026.rds")
# df_sim <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results_sim_noPseudoCount_06-29-2026.rds")

# +1 for all variables
df_mel <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results_mel_1PseudoAll_06-29-2026.rds")
df_sim <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results_sim_1PseudoAll_06-29-2026.rds")

# +1 for SPns and SPs only
# df_mel <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results_mel_1PseudoSPxOnly_06-29-2026.rds")
# df_sim <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/test/adaptedMK_results_sim_1PseudoSPxOnly_06-29-2026.rds")


df_mel[, background := "mel"]
df_sim[, background := "sim"]

df_mel <- melt(
    df_mel, 
    id.vars = 
    c("gene"
    # , "mel_Ps", "mel_Pns",
    # "tsp_SPs", "tsp_SPns", 
    # "conv_SPs", "conv_SPns", 
    # "all_SPs", "all_SPns"
    ),
    measure.vars = c("alpha_tsp", "alpha_conv", "alpha_all"),
    variable.name = "category",
    value.name = "alpha"
)
df_mel[, category := fifelse(category == "alpha_tsp",  "TSP",
                  fifelse(category == "alpha_conv", "CONV",
                          "ALL"))]


df_sim <- melt(
    df_sim, 
    id.vars = 
    c("gene"
    # , "sim_Ps", "sim_Pns",
    # "tsp_SPs", "tsp_SPns", 
    # "conv_SPs", "conv_SPns", 
    # "all_SPs", "all_SPns"
    ),
    measure.vars = c("alpha_tsp", "alpha_conv", "alpha_all"),
    variable.name = "category",
    value.name = "alpha"
)
df_sim[, category := fifelse(category == "alpha_tsp",  "TSP",
                  fifelse(category == "alpha_conv", "CONV",
                          "ALL"))]



df_mel[, background := "mel"]
df_sim[, background := "sim"]
all_df <- rbindlist(list(df_mel, df_sim))

### box plot of mel vs sim alpha values, reported for all categories (ALL, CONV, TSP)
ggplot(all_df, aes(x = background, y = alpha, fill = background)) + geom_boxplot() +
  facet_wrap(~category) + theme_classic() +
  coord_cartesian(ylim = c(-5, 3))


### scatter plot : species-specific conv vs tsp alpha comparison... is TSP alpha > CONV alpha?
ggplot(gene_wide, aes(x = sim_CONV, y = sim_TSP)) +
  geom_point(alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = "sim: TSP vs CONV α",
       x = "α_CONV",
       y = "α_TSP")

### scatter plot : category-specific mel vs sim alpha comparison... is mel alpha similar to sim alpha? 
ggplot(gene_wide, aes(x = mel_TSP, y = sim_TSP)) +
  geom_point(alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = "TSP α: sim vs mel",
       x = "mel_TSP",
       y = "sim_TSP")



