library(data.table)
library(ggplot2)

get_plot_dt <- function(dt){
    dt[, group := 
    fifelse(
        classification%in%tsp_classes, "TSP Candidates", 
        fifelse(
            classification%in%conv_classes, "Convergent",
            fifelse(
                !is.na(classification) & !classification%in%c(tsp_classes, conv_classes), 
                "Other Shared", NA_real_) 
        )
    )]
    class_counts <- dt[!is.na(group) & group != NA_real_, .N, by = group]
    table(dt$group, useNA = "ifany")

    
    ggplot(class_counts, aes(x = group, y = N)) +
        geom_col() +
        geom_text(
            aes(
                label = N,
                y = pmin(N, 70000)
            ),
            angle = 45,
            hjust = 0,
            vjust = -1, 
            size = 2
        ) +
        coord_cartesian(ylim = c(0, 80000)) +
        labs(
            x = "Classification",
            y = "Number of SNPs",
            title = "Distribution of SNPs by classification"
        ) +
        theme_classic()

}

### check ns and syn shared groupings **** 

cand_classes <- c("A", "B", "F", "G", "O", "P", "X", "Y")

tsp_classes <- c("A", "B")
conv_classes <- c("F", "G", "O", "P", "X", "Y")

NS_tsp_classes <- c("A")
NS_conv_classes <- c("F", "O", "X")
NS_shared_other <- c("C", "H", "K", "Q", "T")

syn_tsp_classes <- c("B")
syn_conv_classes <- c("G", "P", "Y")
syn_shared_other <- c("L", "U")

other_shared_other <- c("D", "E", "I", "J", "M", "N", "R", "S", "V", "W")

poolseq_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/poolseq_dest/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.shared.bothMelSim.classed.rds")
inbred_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.11_2026.shared.bothMelSim.0.05MAFfilter.classed.rds")

poolseq_plot_dt <- copy(poolseq_dt)
inbred_plot_dt  <- copy(inbred_dt)

# Add dataset labels
poolseq_plot_dt[, dataset := "PoolSeq"]
inbred_plot_dt[, dataset := "Inbred"]

plot_dt <- rbindlist(list(inbred_plot_dt, poolseq_plot_dt),
    use.names = TRUE, fill = TRUE)

plot_dt[, group := fcase(
    classification %in% tsp_classes,  "TSP Candidates",
    classification %in% conv_classes, "Convergent",
    !is.na(classification),          "Other Shared",
    default = NA_character_
)]

# Assign Syn vs NS
plot_dt[, effect := fcase(
    classification %in% NS_tsp_classes |
        classification %in% NS_conv_classes |
        classification %in% NS_shared_other, "NS",
    classification %in% syn_tsp_classes |
        classification %in% syn_conv_classes |
        classification %in% syn_shared_other, "Syn",
    default = "Other"
)]

# set ordering
plot_dt[, dataset := factor(dataset, levels = c("Inbred", "PoolSeq"))]
plot_dt[, group := factor(group, levels = c("TSP Candidates", "Convergent"
# , "Other Shared"
))]
plot_dt[, effect := factor(effect, levels = c("Syn", "NS", "Other"))]

plot_counts <- plot_dt[!is.na(group), .N, by = .(dataset, group, effect)]

ggplot(plot_counts, aes(x = group, y = N, fill = effect)) +
    geom_col(width = 0.7) +
    facet_grid(. ~ dataset, scales = "fixed", switch = "x") +
    scale_fill_manual(
        values = c(
            "Syn" = "#9ad9db",
            "NS" = "#f49892",
            "Other" = "grey75"
        ), drop = FALSE) +
    labs(
        x = NULL,
        y = "Number of SNPs",
        fill = "Effect"
        # ,
        # title = "Composition of Shared Polymorphism Classifications"
    ) +
    theme_classic() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12),
        axis.text.x = element_text(
            angle = 30,
            hjust = 1
        ),
        legend.position = "right"
    )




##################################################################################
### OLD 






# load rds table final_dt

final_dt[, class_plot := fifelse(
    !is.na(classification),
    classification,
    fifelse(
        !is.na(variant.id_mel) & is.na(variant.id_sim),
        "mel_only",
        fifelse(
            is.na(variant.id_mel) & !is.na(variant.id_sim),
            "sim_only",
            NA_character_
        )
    )
)]

table(final_dt$class_plot, useNA = "ifany")

library(ggplot2)

class_counts <- final_dt[
    !is.na(class_plot) & classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y"),
    .N,
    by = class_plot
]

class_counts[, class_plot := factor(
    class_plot,
    levels = c(LETTERS[1:25], "mel_only", "sim_only")
)]

ggplot(class_counts, aes(x = class_plot, y = N)) +
    geom_col() +
    geom_text(
        aes(
            label = N,
            y = pmin(N, 70000)
        ),
        angle = 45,
        hjust = 0,
        vjust = -1, 
        size = 2
    ) +
    coord_cartesian(ylim = c(0, 80000)) +
    labs(
        x = "Classification",
        y = "Number of SNPs",
        title = "Distribution of SNPs by classification"
    ) +
    theme_classic()