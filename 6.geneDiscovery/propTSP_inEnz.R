library(ggplot2)
library(data.table)

cand_classes <- c("A", "B", "F", "G", "O", "P", "X", "Y")
tsp_classes <- c("A", "B")
conv_classes <- c("F", "G", "O", "P", "X", "Y")

enz <- fread("/scratch/ejy4bu/drosophila/Dmel_enzyme_data_fb_2026_03.tsv", skip=4)
enz.ag <- enz[,list(enz=.N),gene_id]
enzyme_genes <- unique(enz$gene_id)

# load tables (including non-shared and shared variation, both TSP/conv and different)
poolseq_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/poolseq_dest/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.shared.bothMelSim.classed.rds")
inbred_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.11_2026.shared.bothMelSim.0.05MAFfilter.classed.rds")
gene_col <- "gene_id_fbgn_fromGFF"

shared_dt <- merge(
    inbred_dt[, .(chr, pos, classification, gene_id_fbgn_fromGFF)],
    poolseq_dt[, .(chr, pos, classification, gene_id_fbgn_fromGFF)],
    by = c("chr", "pos", "classification"),
    suffixes = c("_inbred", "_poolseq")
)
nrow(shared_dt[gene_id_fbgn_fromGFF_inbred != gene_id_fbgn_fromGFF_poolseq])
shared_dt[, gene_id_fbgn_fromGFF := gene_id_fbgn_fromGFF_inbred]
shared_dt <- unique(shared_dt[!is.na(gene_id_fbgn_fromGFF), .(chr, pos, classification, gene_id_fbgn_fromGFF)])
shared_dt[, is_enzyme := gene_id_fbgn_fromGFF %in% enzyme_genes]


inbred_tsp <- inbred_dt[classification %in% tsp_classes]
poolseq_tsp <- poolseq_dt[classification %in% tsp_classes]

shared_tsp <- merge(
    inbred_tsp[, .(chr, pos, classification, gene_id_fbgn_fromGFF)],
    poolseq_tsp[, .(chr, pos, classification, gene_id_fbgn_fromGFF)],
    by = c("chr", "pos", "classification"),
    suffixes = c("_inbred", "_poolseq")
)

nrow(shared_tsp[gene_id_fbgn_fromGFF_inbred != gene_id_fbgn_fromGFF_poolseq])
shared_tsp[, gene_id_fbgn_fromGFF := gene_id_fbgn_fromGFF_inbred]
shared_tsp <- unique(shared_tsp[!is.na(gene_id_fbgn_fromGFF), .(chr, pos, classification, gene_id_fbgn_fromGFF)])

# Enzyme annotation
inbred_dt[, is_enzyme := gene_id_fbgn_fromGFF %in% enzyme_genes]
poolseq_dt[, is_enzyme := gene_id_fbgn_fromGFF %in% enzyme_genes]
shared_tsp[, is_enzyme := gene_id_fbgn_fromGFF %in% enzyme_genes]

enzyme_analysis <- function(dt, tsp_classes, group_name) {
    tsp <- dt[classification %in% tsp_classes]
    non_tsp <- dt[!classification %in% tsp_classes]

    # 2x2 table: TSP vs non-TSP, enzyme vs non-enzyme
    table <- matrix(
        c(
            sum(tsp$is_enzyme, na.rm = TRUE),
            sum(!tsp$is_enzyme, na.rm = TRUE),
            sum(non_tsp$is_enzyme, na.rm = TRUE),
            sum(!non_tsp$is_enzyme, na.rm = TRUE)
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
            Variation = c("TSP", "NonTSP"),
            Gene = c("Enzyme", "NonEnzyme")
        )
    )

    fisher <- fisher.test(table)

    data.table(
        group = group_name,
        TSP_variants = nrow(tsp),
        TSP_in_enzymes = sum(tsp$is_enzyme, na.rm = TRUE),
        TSP_in_nonenzyme = sum(!tsp$is_enzyme, na.rm = TRUE),
        prop_TSP_in_enzymes = mean(tsp$is_enzyme, na.rm = TRUE),
        nonTSP_variants = nrow(non_tsp),
        nonTSP_in_enzymes = sum(non_tsp$is_enzyme, na.rm = TRUE),
        prop_nonTSP_in_enzymes = mean(non_tsp$is_enzyme, na.rm = TRUE),
        fold_enrichment = mean(tsp$is_enzyme, na.rm = TRUE) /
            mean(non_tsp$is_enzyme, na.rm = TRUE),
        odds_ratio = unname(fisher$estimate),
        p_value = fisher$p.value
    )
}


enzyme_results <- rbindlist(list(

    enzyme_analysis(inbred_dt, tsp_classes, "Inbred"),
    enzyme_analysis(poolseq_dt, tsp_classes, "PoolSeq"),
    enzyme_analysis(shared_dt, tsp_classes, "Shared")
))

enzyme_results

enzyme_results_pretty <- copy(enzyme_results)

enzyme_results_pretty[, `:=`(
    prop_TSP_in_enzymes = round(prop_TSP_in_enzymes * 100, 2),
    prop_nonTSP_in_enzymes = round(prop_nonTSP_in_enzymes * 100, 2),
    fold_enrichment = round(fold_enrichment, 2),
    odds_ratio = round(odds_ratio, 2),
    p_value = signif(p_value, 3)
)]

enzyme_results_pretty

### plot
enzyme_plot <- melt(
    enzyme_results,
    id.vars = "group",
    measure.vars = c(
        "prop_TSP_in_enzymes",
        "prop_nonTSP_in_enzymes"
    ),
    variable.name = "category",
    value.name = "proportion"
)

enzyme_plot[, category := factor(
    category,
    levels = c("prop_TSP_in_enzymes", "prop_nonTSP_in_enzymes"),
    labels = c("TSPs", "Non-TSP variants")
)]

ggplot(enzyme_plot, aes(x = group, y = proportion, fill = category)) +
    geom_col(position = "dodge") +
    scale_y_continuous(labels = scales::percent) +
    labs(
        title = "Proportion of Variants in Enzyme-Encoding Genes",
        x = NULL,
        y = "Proportion of variants in enzyme genes",
        fill = NULL
    ) +
    theme_classic()


shared_dt[, .(
    N = .N,
    enzyme_variants = sum(is_enzyme),
    prop_enzyme = mean(is_enzyme)
), by = classification][order(classification)]

shared_dt[, .(
    N = .N,
    prop_enzyme = mean(is_enzyme)
), by = .(group = fcase(
    classification %in% tsp_classes, "TSP",
    classification %in% conv_classes, "Convergent",
    default = "Other"
))]
#         group      N prop_enzyme
#        <char>  <int>       <num>
# 1:      Other 145708   0.3045612
# 2:        TSP   3001   0.3298900
# 3: Convergent   1180   0.3186441


### old - incorrectly measured prop of genes with TSPs vs genes without TSPs
# Bg genes – all genes in dataset (inbred | poolseq | shared genes)
# Measures prop of enzymes with TSPs vs all enzymes in that dataset 

gene_enzyme_analysis <- function(tsp_dt, background_genes, group_name) {
    tsp_genes <- unique(tsp_dt[!is.na(gene_id_fbgn_fromGFF), gene_id_fbgn_fromGFF])

    background_genes <- unique(background_genes[!is.na(background_genes)])
    tsp_genes <- intersect(tsp_genes, background_genes)

    is_tsp <- background_genes %in% tsp_genes
    is_enzyme <- background_genes %in% enzyme_genes

    # Explicitly ordered: rows = TSP / non-TSP, columns = enzyme / non-enzyme
    table <- matrix(
        c(
            sum(is_tsp & is_enzyme),
            sum(is_tsp & !is_enzyme),
            sum(!is_tsp & is_enzyme),
            sum(!is_tsp & !is_enzyme)
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
            TSP = c("TSP", "NonTSP"),
            Enzyme = c("Enzyme", "NonEnzyme")
        )
    )

    fisher <- fisher.test(table)

    prop_TSP_enz <- mean(is_enzyme[is_tsp])
    prop_bg_enz <- mean(is_enzyme[!is_tsp])

    data.table(
        group = group_name,
        background_genes = length(background_genes),
        TSP_genes = sum(is_tsp),
        TSP_enzyme_genes = sum(is_tsp & is_enzyme),
        TSP_nonenzyme_genes = sum(is_tsp & !is_enzyme),
        prop_TSP_enzymes = prop_TSP_enz,
        prop_background_enzymes = prop_bg_enz,
        fold_enrichment = prop_TSP_enz / prop_bg_enz,
        odds_ratio = unname(fisher$estimate),
        p_value = fisher$p.value
    )
}

# Gene backgrounds
inbred_background <- unique(inbred_dt[!is.na(gene_id_fbgn_fromGFF), gene_id_fbgn_fromGFF])

poolseq_background <- unique(poolseq_dt[!is.na(gene_id_fbgn_fromGFF), gene_id_fbgn_fromGFF])

shared_background <- intersect( inbred_background, poolseq_background)

# Run analyses
enzyme_results <- rbindlist(list(
    gene_enzyme_analysis(
        inbred_tsp,
        inbred_background,
        "Inbred TSPs"
    ),

    gene_enzyme_analysis(
        poolseq_tsp,
        poolseq_background,
        "PoolSeq TSPs"
    ),

    gene_enzyme_analysis(
        shared_tsp,
        shared_background,
        "Shared TSPs"
    )
))

enzyme_results

enzyme_results[, `:=`(
    prop_TSP_enzymes = round(prop_TSP_enzymes * 100, 2),
    prop_background_enzymes = round(prop_background_enzymes * 100, 2),
    fold_enrichment = round(fold_enrichment, 2),
    odds_ratio = round(odds_ratio, 2),
    p_value = signif(p_value, 3)
)]

enzyme_results

### PLOT
enzyme_plot <- melt(
    enzyme_results,
    id.vars = "group",
    measure.vars = c(
        "prop_TSP_enzymes",
        "prop_background_enzymes"
    ),
    variable.name = "category",
    value.name = "proportion"
)

enzyme_plot[, category := factor(
    category,
    levels = c("prop_TSP_enzymes", "prop_background_enzymes"),
    labels = c("TSP-containing genes", "Background genes")
)]

ggplot(enzyme_plot, aes(x = group, y = proportion, fill = category)) +
    geom_col(position = "dodge") +
    labs(
        title = "Enrichment of TSPs in Enzyme-Encoding Genes",
        x = NULL,
        y = "Proportion of genes encoding enzymes",
        fill = NULL
    ) +
    theme_classic()