library(data.table)
library(ggplot2)

geva <- fread("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/AlleleAges.VA.cm_GEVA.txt")


rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total rows in shared table")


# extract chr from id 
geva[, chr := tstrsplit(id, ".", fixed = TRUE)[[1]]]

# table(geva$chr, useNA = "ifany")

setnames(geva, "position", "pos")

# merge geva and rds file 

shared_dt <- merge(
    shared_dt, geva[, .(chr, pos, PostMedian, PostMode)],
    by = c("chr", "pos"),
    all.x=T
)

# number of rows that are missing age data.table
sum(!is.na(shared_dt$PostMedian))

# out_rds <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed_geva.rds")

# saveRDS(shared_dt, out_rds)
# message("classified RDS written to: ", out_rds)

# # CSV SUBSET
# subset_table <- shared_dt[1:500, ]
# # table(subset_table$classification, useNA = "ifany")
# fwrite(subset_table, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed_age_test500.csv")

plot_dt <- shared_dt[!is.na(PostMedian)]

plot <- ggplot(plot_dt, aes(x = classification, y = PostMedian, color=classification)) + 
    geom_point(alpha=0.7) + 
    # geom_boxplot(outlier.alpha = 0.3) +
    theme_classic() + 
    labs(
        x = "Classification", 
        y = "Allele Age (generations)",
        title = "GEVA Allele Age Estimates by Polymorphism Class",
    ) + 
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        text = element_text(size = 15), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
    ) + 
    scale_y_log10() +
    scale_color_manual(
        values = c(
            "A" = "deepskyblue1",
            "B" = "deepskyblue1",
            "C" = "firebrick",
            "D" = "firebrick",
            "E" = "firebrick",

            "F" = "steelblue",
            "G" = "steelblue",

            "H" = "firebrick",
            "I" = "firebrick",
            "J" = "firebrick",
            "K" = "firebrick",
            "L" = "firebrick",
            "M" = "firebrick",
            "N" = "firebrick",

            "O" = "steelblue",
            "P" = "steelblue",

            "Q" = "firebrick",
            "R" = "firebrick",
            "S" = "firebrick",
            "T" = "firebrick",
            "U" = "firebrick",
            "V" = "firebrick",
            "W" = "firebrick",

            "X" = "steelblue",
            "Y" = "steelblue"

            # "Convergent" = "steelblue",
            # "Divergent" = "firebrick"
        )
    ) + 
    annotate(
        "label",
        x = "T",
        y = max(plot_dt$PostMode) * 0.95,
        label = "Blue = Convergent\nRed = Divergent",
        fill = "white",
        color = "black",
        size = 4
    )


out_dir <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/figures/alleleAge")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
    filename = file.path(out_dir, paste0("allClasses_median_scatter_logYaxis.png")),
    plot = plot,
    width = 8, height=6,
    dpi=300
)