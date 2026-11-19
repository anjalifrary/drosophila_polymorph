library(data.table)
library(ggplot2)

out_dir <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/figures/classificationPlots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

class_table <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/classification/classification_table_final.csv")
class_table <- read.csv(class_table)

names(class_table)
setnames(class_table,
         "Convergent.Divergent",
         "conv_class")

plot <- ggplot(class_table, aes(x=Classification, y = Count, fill = conv_class)) + 
    geom_bar(stat = "identity", color = "black", width = 0.5) + 
    theme_bw() + 
    labs(
        x = "Classification", 
        y = "Number of variants",
        title = "Distribution of SNPs by class",
        fill = "Class"
    ) + 
    theme_classic() + 
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        text = element_text(size = 15), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
    ) + 
    scale_fill_manual(
        values = c(
            "Convergent" = "steelblue",
            "Divergent" = "firebrick"
        )
    ) + 
    annotate(
        "label",
        x = 22,
        y = max(classes$Count) * 0.95,
        label = "Blue = Convergent\nRed = Divergent",
        fill = "white",
        color = "black",
        size = 4
    )
    
# ggsave(
#     filename = file.path(out_dir, paste0("allClasses_barPlot_groupedByConvergence.png")),
#     plot = plot,
#     width = 8, height=6,
#     dpi=300
# )


### new 11/19/2026, ran after propTSP_inEnz.R (using those objects)

# Count variants by classification group
count_classes <- function(dt, dataset) {
    dt[, group := fcase(
        classification %in% tsp_classes, "TSP (A/B)",
        classification %in% conv_classes, "Convergent",
        default = "Other"
    )]

    dt[, .(Count = .N), by = group][, dataset := dataset]
}

class_counts <- rbindlist(list(
    count_classes(copy(inbred_dt), "Inbred"),
    count_classes(copy(poolseq_dt), "PoolSeq"),
    count_classes(copy(shared_dt), "Shared")
))

class_counts[, group := factor(
    group,
    levels = c("TSP (A/B)", "Convergent", "Other")
)]

class_counts[, dataset := factor(
    dataset,
    levels = c("Inbred", "PoolSeq", "Shared")
)]

class_counts[, prop := Count / sum(Count), by = dataset]

# Plot
# class_plot <- ggplot(class_counts, aes(x = dataset, y = Count, fill = group)) +
class_plot <- ggplot(class_counts, aes(x = dataset, y = prop, fill = group)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c(
        "TSP (A/B)" = "brown1",
        "Convergent" = "deepskyblue1",
        "Other" = "grey70"
    )) +
    theme_classic() +
    labs(
        x = NULL,
        # y = "Number of variants",
        y = "Proportion of variants",
        title = "Distribution of SNPs by Classification Group",
        fill = "Classification"
    ) +
    theme(
        panel.grid = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
    )

class_plot

class_plot
