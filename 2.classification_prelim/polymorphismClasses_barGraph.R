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
    
ggsave(
    filename = file.path(out_dir, paste0("allClasses_barPlot_groupedByConvergence.png")),
    plot = plot,
    width = 8, height=6,
    dpi=300
)

# plot_dt[, prop := N / sum(N)]

# plot <- ggplot(plot_dt, aes(x = class_simple, y = N, fill = class_simple)) +
#     geom_bar(stat = "identity", color = "black", width = 0.7) +
#     # scale_y_log10() +
#     theme_bw() + 
#     labs(
#         x = "Same-site classification",
#         y = "Number of variants",
#         title = "Distribution of same-site SNP classes",
#         fill = "Class"
#     ) +
#     theme_classic() + 
#     theme(
#         panel.grid = element_blank(),
#         legend.position = "none",
#         text = element_text(size=16),
#         axis.text.x = element_text(size=14),
#         axis.text.y = element_text(size=14),
#         plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
#     ) +
#     scale_fill_manual(values = c(
#         "A" = "brown1",  # convergent (strong)
#         "C" = "pink1",  # convergent (lighter)
#         "E" = "lightsteelblue2",  # divergent (strong)
#         "B" = "deepskyblue1",
#         "D" = "skyblue3",
#         "F" = "cyan1" 
#     )) 

# ggsave(
#     filename = file.path(out_dir, paste0("sameSiteClassification_bar_pretty.png")),
#     plot = plot,
#     width = 8, height=6,
#     dpi=300
# )

