# https://cran.r-project.org/web/packages/grantham/readme/README.html

# The Grantham distance attempts to provide a proxy for the evolutionary distance 
    # between two amino acids based on three key side chain chemical properties: 
    # composition, polarity and molecular volume. 
# In turn, evolutionary distance is used as a proxy for the impact of missense substitutions. 
# The higher the distance, the more deleterious the substitution is expected to be.


library(grantham)
library(ggplot2)
library(data.table)

grantham_distance(x = "Ser", y = "Phe")
grantham_distance(x = c("Ser", "Arg"), y = c("Phe", "Leu"))

# generates 20x20 AA pairs
aa_pairs <- amino_acid_pairs()
 
# removes Ser-Ser and Ser-Arg Arg-Ser are the same
# aa_pairs <- amino_acid_pairs(keep_self = FALSE, keep_reverses = FALSE)

results <- grantham_distance(x = aa_pairs$x, y = aa_pairs$y)
# results format
#> # A tibble: 400 × 3
#>    x     y         d
#>    <chr> <chr> <dbl>
#>  1 Ser   Ser       0
#>  2 Ser   Arg     110
#>  3 Ser   Leu     145
#>  4 Ser   Pro      74
setnames(results, "d", "grantham_score")


### load data 
snps_0.01 <- readRDS("/scratch/ejy4bu/drosophila/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.candidatesABFGOPXY.0.01MAFfilter.classed.rds")
snps_0.05 <- readRDS("/scratch/ejy4bu/drosophila/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.candidatesABFGOPXY.0.05MAFfilter.classed.rds")

snps_0.01 <- readRDS("/scratch/ejy4bu/drosophila/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.shared.bothMelSim.0.01MAFfilter.classed.rds")
snps_0.05 <- readRDS("/scratch/ejy4bu/drosophila/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.shared.bothMelSim.0.05MAFfilter.classed.rds")

tsps <- snps_0.05[classification%in%c("A", "B")]
genes <- unique(tsps[, gene_id_fbgn_fromGFF])
gene_snps <- tsps[
    ,
    .(
        total_snps = .N,
        A_snps = sum(classification == "A"),
        B_snps = sum(classification == "B")
    ),
    by = gene_id_fbgn_fromGFF
]

setorder(gene_snps, -total_snps)
gene_snps[, .(
     total = sum(total_snps),
     A = sum(A_snps),
     B = sum(B_snps)
)]

enz <- fread("/scratch/ejy4bu/drosophila/Dmel_enzyme_data_fb_2026_03.tsv", skip=4)
enz.ag <- enz[,list(enz=.N),gene_id]
tsps <- merge(tsps, enz.ag, by.x="gene_id_fbgn_fromGFF", by.y="gene_id", all.x=T)
enzymes <- tsps[enz > 0]
enzymes <- merge(enzymes, gene_snps, by = "gene_id_fbgn_fromGFF", all.x=T)
candidate_ns_enzymes <- enzymes[A_snps>0]
candidate_ns_enzymes <- candidate_ns_enzymes[af_mel>0.25 & af_sim>0.25]
candidate_ns_snps <- candidate_ns_enzymes[classification=="A"]

#### add grantham scores
snpDT <- snps_0.05

snpDT <- merge(snpDT, results, by.x=c("aa_ref_mel", "aa_alt_mel"), by.y=c("x", "y"), all.x=T)
setnames(snpDT, "grantham_score", "grantham_score_mel")
snpDT <- merge(snpDT, results, by.x=c("aa_ref_sim", "aa_alt_sim"), by.y=c("x", "y"), all.x=T)
setnames(snpDT, "grantham_score", "grantham_score_sim")

snpDT[, grantham_score := fifelse(!is.na(grantham_score_mel), grantham_score_mel, grantham_score_sim)]
snpDT[, grantham_score_mel := NULL]
snpDT[, grantham_score_sim := NULL]


hist(snpDT[classification=="A", ]$grantham_score, breaks=100, main="Distribution of Grantham Scores - Class A", xlab = "Grantham score")
candidate_ns_snps <- merge(candidate_ns_snps, snpDT[, .(chr, pos, grantham_score)], by=c("chr", "pos"), all.x=T)
hist(candidate_ns_snps$grantham_score, breaks=100, main="Distribution of Grantham Scores - Class A enzyme candidates", xlab = "Grantham score")


### define comparison classes
grantham_A <- snpDT[classification == "A", .(grantham_score, class_group = "Class A")]
grantham_tsp_conv <- snpDT[classification %in% c("B", "F", "G", "O", "P", "X", "Y"), 
    .(grantham_score, class_group = "TSP & Conv")
]

grantham_all_classed <- snpDT[
    !is.na(classification),
    .(grantham_score, class_group = "All classed")
]

grantham_mel_only <- snpDT[
    !is.na(variant.id_mel) & is.na(variant.id_sim),
    .(grantham_score, class_group = "Mel-only")
]
grantham_sim_only <- snpDT[
    is.na(variant.id_mel) & !is.na(variant.id_sim),
    .(grantham_score, class_group = "Sim-only")
]

grantham_plot <- rbindlist(list(
    grantham_A,
    grantham_tsp_conv,
    grantham_all_classed,
    grantham_mel_only,
    grantham_sim_only
))

grantham_plot <- grantham_plot[
    !is.na(grantham_score) &
    # is.finite(grantham_score) &
    grantham_score > 0
]

grantham_plot[, class_group := factor(
    class_group,
    levels = c(
        "Class A",
        "TSP & Conv",
        "All classed",
        "Mel-only",
        "Sim-only"
    )
)]

### boxplot
ggplot(
    grantham_plot,
    aes(x = class_group, y = grantham_score)
) +
    geom_boxplot(outlier.alpha = 0.1) +
    labs(
        title = "Grantham scores by classification",
        x = "Classification",
        y = "Grantham score"
    ) +
    theme_classic()

### hist, by group
ggplot(
    grantham_plot,
    aes(x = grantham_score)
) +
    geom_histogram(bins = 50) +
    facet_wrap(~class_group, nrow = 1) +
    labs(
        title = "Distribution of Grantham scores",
        x = "Grantham score",
        y = "Number of SNPs"
    ) +
    theme_classic()

### mean
grantham_summary <- grantham_plot[
    ,
    .(
        mean_grantham = mean(grantham_score),
        sd_grantham = sd(grantham_score),
        N = .N
    ),
    by = class_group
]

### mean
ggplot(
    grantham_summary,
    aes(x = class_group, y = mean_grantham)
) +
    geom_col(fill = "steelblue") +
    labs(
        title = "Mean Grantham score by classification",
        x = "Classification",
        y = "Mean Grantham score"
    ) +
    theme_classic()

### median
grantham_summary <- grantham_plot[
    ,
    .(
        mean_grantham = mean(grantham_score),
        median_grantham = median(grantham_score),
        N = .N
    ),
    by = class_group
]

ggplot(
    grantham_summary,
    aes(x = class_group, y = median_grantham)
) +
    geom_col(fill = "steelblue") +
    labs(
        title = "Median Grantham score by classification",
        x = "Classification",
        y = "Median Grantham score"
    ) +
    theme_classic()







    

ggplot(
    grantham_plot,
    aes(x = class_group, y = grantham_score)
) +
    geom_boxplot(outlier.alpha = 0.1) +
    labs(
        title = "Grantham scores by classification",
        x = "Classification group",
        y = "Grantham score"
    ) +
    theme_classic()