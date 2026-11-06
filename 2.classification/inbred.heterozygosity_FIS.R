library(data.table)

###############################################################################################
### Heterozygosity & FIS
###############################################################################################

mel_snpdata <- readRDS("/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/mel_site_RD.rds")
sim_snpdata <- readRDS("/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/sim_site_RD.rds")
shared <- readRDS("/scratch/ejy4bu/drosophila/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.shared.0.05MAFfilter.classed.rds")

shared_het <- merge(shared, mel_snpdata[, .(chr, pos, prop_het, fis)], by=c("chr", "pos"), all.x=T)
setnames(shared_het, c("prop_het", "fis"), c("prop_het_mel", "fis_mel"))
shared_het <- merge(shared_het, sim_snpdata[, .(chr, pos, prop_het, fis)], by=c("chr", "pos"), all.x=T)
setnames(shared_het, c("prop_het", "fis"), c("prop_het_sim", "fis_sim"))

shared_het[, .(
    total = .N,
    mel_fis_present = sum(!is.na(fis_mel)),
    mel_fis_missing = sum(is.na(fis_mel)),
    sim_fis_present = sum(!is.na(fis_sim)),
    sim_fis_missing = sum(is.na(fis_sim))
)]

# Compare variant IDs against FIS availability
shared_het[, .(
    mel_variant_but_missing_fis =
        sum(!is.na(variant.id_mel) & is.na(fis_mel)),
    sim_variant_but_missing_fis =
        sum(!is.na(variant.id_sim) & is.na(fis_sim))
)]

library(data.table)
library(ggplot2)

# Define classification groups
shared_het_plot <- shared_het[!is.na(classification)]
shared_het_plot[, class_group := fifelse(
    classification %in% c("A", "B"),
    "tsp",
    "Not tsp"
)]

fis_long <- melt(
    shared_het_plot,
    id.vars = "class_group",
    measure.vars = c("fis_mel", "fis_sim"),
    variable.name = "species",
    value.name = "FIS"
)
fis_long[, species := fifelse(
    species == "fis_mel",
    "D. melanogaster",
    "D. simulans"
)]
fis_long <- fis_long[
    !is.na(FIS) & is.finite(FIS)
]

ggplot(fis_long, aes(x = class_group, y = FIS)) +
    geom_boxplot(outlier.alpha = 0.1) +
    facet_wrap(~species, nrow = 1) +
    labs(
        title = expression("Site-level " * F[IS] * " by classification"),
        x = "Classification group",
        y = expression(F[IS])
    ) +
    theme_classic()


fis_prop_negative <- melt(
    shared_het_plot,
    id.vars = "class_group",
    measure.vars = c("fis_mel", "fis_sim"),
    variable.name = "species",
    value.name = "FIS"
)[
    !is.na(FIS) & is.finite(FIS),
    .(
        n_sites = .N,
        n_negative = sum(FIS < 0),
        prop_negative = mean(FIS < 0),
        percent_negative = 100 * mean(FIS < 0)
    ),
    by = .(species, class_group)
][order(species, class_group)]

library(ggplot2)

ggplot(
    fis_prop_negative,
    aes(x = species, y = prop_negative, fill = class_group)
) +
    geom_col(position = "dodge") +
    labs(
        x = "Classification",
        y = "Proportion with FIS < 0",
        fill = "Species"
    ) +
    theme_classic()


    prop_het_long <- melt(
    shared_het_plot,
    id.vars = "class_group",
    measure.vars = c("prop_het_mel", "prop_het_sim"),
    variable.name = "species",
    value.name = "prop_het"
)

prop_het_long[, species := fifelse(
    species == "prop_het_mel",
    "D. melanogaster",
    "D. simulans"
)]

prop_het_long <- prop_het_long[
    !is.na(prop_het) & is.finite(prop_het)
]

ggplot(prop_het_long, aes(x = class_group, y = prop_het)) +
    geom_boxplot(outlier.alpha = 0.1) +
    facet_wrap(~species, nrow = 1) +
    labs(
        title = "Observed heterozygosity by classification",
        x = "Classification group",
        y = "Proportion heterozygous"
    ) +
    theme_classic()

prop_het_summary <- prop_het_long[
    ,
    .(
        n_sites = .N,
        mean_prop_het = mean(prop_het),
        sd_prop_het = sd(prop_het)
    ),
    by = .(species, class_group)
]

ggplot(
    prop_het_summary,
    aes(x = species, y = mean_prop_het, fill = class_group)
) +
    geom_col(position = "dodge") +
    labs(
        title = "Mean observed heterozygosity by classification",
        x = "Species",
        y = "Mean proportion heterozygous",
        fill = "Classification"
    ) +
    theme_classic()

