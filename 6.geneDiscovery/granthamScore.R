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

grantham_lookup  <- grantham_distance(x = aa_pairs$x, y = aa_pairs$y)
# results format
#> # A tibble: 400 × 3
#>    x     y         d
#>    <chr> <chr> <dbl>
#>  1 Ser   Ser       0
#>  2 Ser   Arg     110
#>  3 Ser   Leu     145
#>  4 Ser   Pro      74
setnames(grantham_lookup , "d", "grantham_score")
setDT(grantham_lookup)
grantham_lookup[, pair := paste(x, y, sep = "_")]

grantham_lookup <- grantham_lookup[, .(pair, grantham_score)]


# Function to add Grantham scores to a SNP table
add_grantham <- function(dt) {
    dt <- copy(dt)

    # Mel substitution score
    mel_pair <- paste(dt$aa_ref_mel, dt$aa_alt_mel,sep = "_")

    # Sim substitution score
    sim_pair <- paste(dt$aa_ref_sim,dt$aa_alt_sim,sep = "_")

    score_lookup <- setNames(grantham_lookup$grantham_score,grantham_lookup$pair)

    dt[, grantham_score_mel := unname(score_lookup[mel_pair])]
    dt[, grantham_score_sim := unname(score_lookup[sim_pair])]
    # stopifnot(all(dt$grantham_score_mel == dt$grantham_score_sim, na.rm = TRUE))

    dt[, grantham_score := grantham_score_mel]

    # dt[, c("grantham_score_mel", "grantham_score_sim") := NULL]

    dt
}

gene_level_analysis <- function(dt, dataset_name) {

    gene_TSP <- dt[
        !is.na(gene_id_fbgn_fromGFF),
        .(
            total_TSP = sum(classification %in% c("A", "B")),
            A_snps = sum(classification == "A"),
            B_snps = sum(classification == "B")
        ),
        by = gene_id_fbgn_fromGFF
    ]

    gene_grantham <- dt[
        classification == "A" &
        !is.na(grantham_score) &
        grantham_score > 0 &
        !is.na(gene_id_fbgn_fromGFF),
        .(
            n_A_scored = .N,
            mean_grantham_A = mean(grantham_score),
            max_grantham_A = max(grantham_score)
        ),
        by = gene_id_fbgn_fromGFF
    ]

    gene_metrics <- merge(
        gene_TSP,
        gene_grantham,
        by = "gene_id_fbgn_fromGFF",
        all.x = TRUE
    )

    gene_metrics[
        is.na(n_A_scored),
        n_A_scored := 0L
    ]

    stopifnot(all(dt$A_snps == dt$n_A_scored, na.rm = TRUE))
    gene_metrics[, n_A_scored := NULL]

    gene_metrics[, dataset := dataset_name]

    setorder(gene_metrics, -total_TSP)

    gene_metrics
}


# load tables (including non-shared and shared variation, both TSP/conv and different)
poolseq_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/poolseq_dest/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.shared.bothMelSim.classed.rds")
inbred_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.11_2026.shared.bothMelSim.0.05MAFfilter.classed.rds")
gene_col <- "gene_id_fbgn_fromGFF"

cand_classes <- c("A", "B", "F", "G", "O", "P", "X", "Y")
tsp_classes <- c("A", "B")
conv_classes <- c("F", "G", "O", "P", "X", "Y")

poolseq_dt <- add_grantham(poolseq_dt)
inbred_dt <- add_grantham(inbred_dt)

pooled_tsps <- poolseq_dt[classification%in%tsp_classes,]
inbred_tsps <- inbred_dt[classification%in%tsp_classes,]

poolseq_gene_metrics <- gene_level_analysis(pooled_tsps, "PoolSeq")
inbred_gene_metrics  <- gene_level_analysis(inbred_tsps, "Inbred")

gene_metrics_all <- rbindlist(list(poolseq_gene_metrics, inbred_gene_metrics))

enz <- fread("/scratch/ejy4bu/drosophila/Dmel_enzyme_data_fb_2026_03.tsv", skip=4)
enz.ag <- enz[,list(enz=.N),gene_id]
gene_metrics_all <- merge(gene_metrics_all, enz.ag, by.x="gene_id_fbgn_fromGFF", by.y="gene_id", all.x=T)

make_grantham_plot_dt <- function(dt, dataset_name) {

    rbindlist(list(
        dt[classification == "A",
           .(grantham_score, class_group = "Class A")],

        dt[classification == "B",
           .(grantham_score, class_group = "Class B")],

        dt[classification %in% conv_classes,
           .(grantham_score, class_group = "Convergent")],

        dt[classification %in% c("B", conv_classes),
           .(grantham_score, class_group = "TSP & Conv")],

        dt[!is.na(classification),
           .(grantham_score, class_group = "All classed")]
    ))[, dataset := dataset_name]
}

grantham_plot <- rbindlist(list(
    make_grantham_plot_dt(poolseq_dt, "PoolSeq"),
    make_grantham_plot_dt(inbred_dt, "Inbred")
))

grantham_plot <- grantham_plot[
    !is.na(grantham_score) & grantham_score > 0
]

grantham_plot[, class_group := factor(
    class_group,
    levels = c(
        "Class A",
        "Class B",
        "Convergent",
        "TSP & Conv",
        "All classed"
    )
)]

### boxplots
ggplot(grantham_plot, aes(x = class_group, y = grantham_score)) +
    geom_boxplot(outlier.alpha = 0.1) +
    facet_wrap(~dataset) +
    labs(
        title = "Grantham Scores by Classification",
        x = "Classification",
        y = "Grantham score"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Hists
ggplot(grantham_plot, aes(x = grantham_score)) +
    geom_histogram(bins = 50) +
    facet_grid(dataset ~ class_group) +
    labs(
        title = "Distribution of Grantham Scores",
        x = "Grantham score",
        y = "Number of SNPs"
    ) +
    theme_classic()

# summary stats
grantham_summary <- grantham_plot[
    ,
    .(
        mean_grantham = mean(grantham_score),
        median_grantham = median(grantham_score),
        sd_grantham = sd(grantham_score),
        N = .N
    ),
    by = .(dataset, class_group)
]

grantham_summary

# mean score
ggplot(grantham_summary, aes(x = class_group, y = mean_grantham)) +
    geom_col(fill = "steelblue") +
    facet_wrap(~dataset) +
    labs(
        title = "Mean Grantham Score by Classification",
        x = "Classification",
        y = "Mean Grantham score"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# median
ggplot(grantham_summary, aes(x = class_group, y = median_grantham)) +
    geom_col(fill = "steelblue") +
    facet_wrap(~dataset) +
    labs(
        title = "Median Grantham Score by Classification",
        x = "Classification",
        y = "Median Grantham score"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

### significance tests

grantham_tests <- grantham_plot[
    ,
    {
        A <- grantham_score[class_group == "Class A"]
        Conv <- grantham_score[class_group == "Convergent"]
        TSP_Conv <- grantham_score[class_group == "TSP & Conv"]
        All <- grantham_score[class_group == "All classed"]

        list(
            mean_A = mean(A),
            mean_Conv = mean(Conv),
            mean_TSP_Conv = mean(TSP_Conv),
            mean_All = mean(All),

            difference_A_Conv = mean(A) - mean(Conv),
            p_A_Conv = t.test(A, Conv)$p.value,

            difference_TSP_Conv = mean(TSP_Conv) - mean(Conv),
            p_TSP_Conv = t.test(TSP_Conv, Conv)$p.value,

            difference_All_Conv = mean(All) - mean(Conv),
            p_All_Conv = t.test(All, Conv)$p.value
        )
    },
    by = dataset
]

grantham_tests
### summary table

grantham_summary_table <- rbindlist(list(
    grantham_tests[, .(
        dataset,
        comparison = "Class A vs Convergent",
        mean_group = mean_A,
        mean_convergent = mean_Conv,
        difference = difference_A_Conv,
        p_value = p_A_Conv
    )],

    grantham_tests[, .(
        dataset,
        comparison = "TSP & Conv vs Convergent",
        mean_group = mean_TSP_Conv,
        mean_convergent = mean_Conv,
        difference = difference_TSP_Conv,
        p_value = p_TSP_Conv
    )],

    grantham_tests[, .(
        dataset,
        comparison = "All classed vs Convergent",
        mean_group = mean_All,
        mean_convergent = mean_Conv,
        difference = difference_All_Conv,
        p_value = p_All_Conv
    )]
))

grantham_summary_table[, `:=`(
    mean_group = round(mean_group, 2),
    mean_convergent = round(mean_convergent, 2),
    difference = round(difference, 2),
    p_value = signif(p_value, 3)
)]

grantham_summary_table

### plos difference in mean 

ggplot(grantham_summary_table, aes(x = difference, y = comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 3) +
    facet_wrap(~dataset) +
    labs(
        title = "Differences in Mean Grantham Scores",
        x = "Difference in mean Grantham score (Group - Convergent)",
        y = NULL
    ) +
    theme_classic()


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