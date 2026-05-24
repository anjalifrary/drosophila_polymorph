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
    # scale_y_log10() +
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

### ANOVA

# # make classification categorical
# anova_dt <- shared_dt[!is.na(PostMedian)]
# anova_dt[, classification := factor(classification)]

# # run ANOVA on raw allele ages
# anova_model <- aov(PostMedian ~ classification, data = anova_dt)

# # ANOVA table
# summary(anova_model)


# model <- lm(PostMedian ~ classification, data = anova_dt)
# anova(model)


### from claude
# # Create a binary grouping: A/B vs everything else
# anova_dt[, group := ifelse(classification %in% c("A", "B"), "AB", "Other")]
# anova_dt[, group := factor(group, levels = c("Other", "AB"))]  # Other = reference

# # Simple t-test (or Wilcoxon if non-normal — allele ages often are)
# t.test(PostMedian ~ group, data = anova_dt, alternative = "less")
# # alternative = "less" because Other < AB (you expect AB to be older)

# # Non-parametric version (often better for skewed allele age distributions)
# wilcox.test(PostMedian ~ group, data = anova_dt, alternative = "less")

# # Cohen's d (base R only)
# ab_vals   <- anova_dt[group == "AB",    PostMedian]
# other_vals <- anova_dt[group == "Other", PostMedian]

# pooled_sd <- sqrt((sd(ab_vals)^2 + sd(other_vals)^2) / 2)
# cohens_d  <- (mean(ab_vals) - mean(other_vals)) / pooled_sd
# message("Cohen's d: ", round(cohens_d, 3))

# # Rank-biserial correlation (effect size for Wilcoxon, also base R)
# n_ab    <- length(ab_vals)
# n_other <- length(other_vals)
# W <- wilcox.test(PostMedian ~ group, data = anova_dt)$statistic
# r_rb <- 1 - (2 * W) / (n_ab * n_other)
# message("Rank-biserial r: ", round(r_rb, 3))





# merge tables from all variants (inc not shared) and variants of interest
all_variants <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/merged_tables/quality/all_quality_variants_clean.rds")
classed_variants <- readRDS(rds_file)
all_variants[, classification := NULL]

merged_dt <- merge(
    all_variants, classed_variants,
    by = c("chr", "pos", "ref_mel", "alt_mel",
    "ref_sim", "alt_sim", "nt_change_mel", "nt_change_sim",
    "codon_ref_mel", "codon_alt_mel", "codon_ref_sim", "codon_alt_sim",
    "aa_ref_mel", "aa_alt_mel", "aa_ref_sim", "aa_alt_sim", 
    "af_mel", "af_sim", "effect_mel", "effect_sim", 
    "gene_mel", "gene_sim", "gene_id_mel", "gene_id_sim",
    "aa_pos_mel", "aa_pos_sim"),
    all.x = T
)
message(nrow(merged_dt), " rows in merged table")
message(nrow(all_variants), " rows in all variants table")

# extract chr from id 
geva[, chr := tstrsplit(id, ".", fixed = TRUE)[[1]]]

# table(geva$chr, useNA = "ifany")

setnames(geva, "position", "pos")

merged_dt <- merge(
    merged_dt, geva[, .(chr, pos, PostMedian, PostMode)],
    by = c("chr", "pos"),
    all.x=T
)

# number of rows that are missing age data.table
sum(!is.na(merged_dt$PostMedian))

merged_dt[(classification=="A" | classification=="B"), 
group := "TSP_or_conv"]

merged_dt[(!is.na(ref_mel) & is.na(ref_sim) & is.na(classification)),
group := "mel_only"]

merged_dt[(!is.na(ref_sim) & is.na(ref_mel) & is.na(classification)),
group := "sim_only"]

merged_dt[(classification=="F" | classification=="G" |
    classification=="O" | classification=="P" |
    classification=="X" | classification=="Y"),
    group:="conv_ind"]

merged_dt[(classification=="C" | classification=="D" |
    classification=="E" | classification=="H" |
    classification=="I" | classification=="J" | 
    classification=="K" | classification=="L" |
    classification=="M" | classification=="N" |
    classification=="Q" | classification=="R" |
    classification=="S" | classification=="T" |
    classification=="U" | classification=="V" |
    classification=="W" ),
    group:="div_ind"]

merged_dt[(is.na(classification) & is.na(group)),
    group := "other"]

### go back and select for max 1 polymorphism per codon

table(merged_dt$group, useNA="ifany")


plot_dt <- merged_dt[!is.na(PostMedian)]
plot_dt$group <- factor(plot_dt$group, levels = c("TSP_or_conv", "conv_ind", "div_ind", "mel_only","sim_only", "other"))

plot <- ggplot(plot_dt, aes(x = group, y = PostMedian, color=group)) + 
    # geom_point(alpha=0.7) + 
    geom_boxplot(outlier.alpha = 0.3) +
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
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
    ) + 
    # scale_y_log10() +
    scale_color_manual(
        values = c(
            "TSP_or_conv" = "deepskyblue1",
            "conv_ind" = "steelblue",
            "div_ind" = "firebrick"
        )
    ) + 
    annotate(
        "label",
        x = "other",
        y = max(plot_dt$PostMode) * 0.95,
        label = "Blue = Convergent\nRed = Divergent",
        fill = "white",
        color = "black",
        size = 4
    )

out_dir <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/figures/alleleAge")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
    filename = file.path(out_dir, paste0("allvariants_Nshared_box.png")),
    plot = plot,
    width = 8, height=6,
    dpi=300
)