library(data.table)
library(ggplot2)

out_dir <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/figures/classificationPlots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_variantsOfInterest_qualFiltered.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")

table(shared_dt$classification_sameSite, useNA = "ifany")

# get the number of variants for each class, regardless of the different permutations
plot_dt <- shared_dt[
    !is.na(classification_sameSite),
    .(class_simple = substr(classification_sameSite, 1, 1))
]

# get counts per class
plot_dt <- plot_dt[
    class_simple %in% c("A","B","C","D","E","F"),
    .N,
    by = class_simple
]

# order
plot_dt[, class_simple := factor(
    class_simple,
    levels = c("A","B","C","D","E","F")
)]


plot_dt[, prop := N / sum(N)]

plot <- ggplot(plot_dt, aes(x = class_simple, y = N, fill = class_simple)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    # scale_y_log10() +
    theme_bw() + 
    labs(
        x = "Same-site classification",
        y = "Number of variants",
        title = "Distribution of same-site SNP classes",
        fill = "Class"
    ) +
    theme_classic() + 
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        text = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
    ) +
    scale_fill_manual(values = c(
        "A" = "brown1",  # convergent (strong)
        "C" = "pink1",  # convergent (lighter)
        "E" = "lightsteelblue2",  # divergent (strong)
        "B" = "deepskyblue1",
        "D" = "skyblue3",
        "F" = "cyan1" 
    )) 

ggsave(
    filename = file.path(out_dir, paste0("sameSiteClassification_bar_pretty.png")),
    plot = plot,
    width = 8, height=6,
    dpi=300
)

