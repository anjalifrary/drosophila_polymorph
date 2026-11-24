library(data.table)
library(grantham)
library(ggplot2)
library(stringr)


### load functions from granthamScore.R :
# add_grantham()

library(data.table)

get_longest_transcript_lengths <- function(gff_file) {
    # gff_file <- mel_gff
    gff <- fread(
        cmd = paste("grep -v '^#'", shQuote(gff_file)), 
        sep="\t", header=F, 
        col.names=c(
            "chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"
    ))

    # map chr names to dm6 normal names 
    chr_map <- gff[ feature == "region" & grepl("genome=chromosome", attributes),
        .(chr_gff = chr, chr = sub(".*Name=([^;]+);.*", "\\1", attributes))]

    gff[chr_map, chr_dm6 := i.chr, on = .(chr = chr_gff)]
    unique(gff[!is.na(chr_dm6), .(chr, chr_dm6)]) # X, 2L, 2R, 3L, 3R, 4, Y
    cds_gff <- gff[feature == "CDS" & !is.na(chr_dm6), .(chr_dm6, start, end, strand, phase, attributes)]
    cds_gff[, transcript_id := str_extract(attributes, "(?<=Parent=rna-)[^;,]+")]
    cds_gff[, gene_name := str_extract(attributes, "(?<=;gene=)[^;]+")]
    cds_gff[, gene_id_fbgn := str_extract(attributes, "FBgn[0-9]+")]

    transcript_lengths <- cds_gff[,
        .(transcript_length = sum(end - start + 1)),
        by = .(gene_name, transcript_id, gene_id_fbgn)
    ]

    gene_lengths <- transcript_lengths[, .(
        longest_transcript_length = max(transcript_length),
        longest_transcript_id = transcript_id[
            which.max(transcript_length)
        ]), by = .(gene_id_fbgn, gene_name)
    ]
    gene_lengths
}

summarize_grantham <- function(dt, dataset_name, gene_lengths) {
    dt <- copy(dt)

    # Map gene name and longest CDS length from GFF
    dt <- merge(dt, gene_lengths[, .(gene_id_fbgn, gene_name, longest_transcript_length)],
        by.x = "gene_id_fbgn_fromGFF", by.y = "gene_id_fbgn", all.x = TRUE)

    # Keep valid missense SNPs with matched genes and scores
    dt <- dt[!is.na(gene_id_fbgn_fromGFF) & 
        !is.na(longest_transcript_length) & 
        longest_transcript_length > 0
    ]
    mel_pos <- !is.na(dt$grantham_score_mel) & dt$grantham_score_mel > 0
    sim_pos <- !is.na(dt$grantham_score_sim) & dt$grantham_score_sim > 0
    dt[, grantham_score_combined := NA_real_]
    
    # mel missense only
    dt[mel_pos & !sim_pos, grantham_score_combined := grantham_score_mel]
    # sim missense only
    dt[!mel_pos & sim_pos, grantham_score_combined := grantham_score_sim]
    # both syn
    dt[!is.na(grantham_score_mel) & !is.na(grantham_score_sim) & grantham_score_mel == 0 & grantham_score_sim == 0,
        grantham_score_combined := 0]

    # both missense; same amino acid pair
    dt[mel_pos & sim_pos & grantham_score_mel == grantham_score_sim, grantham_score_combined := grantham_score_mel ]
    # both missense; diff amino acid pairs - idk how to record this
    dt[mel_pos & sim_pos & grantham_score_mel != grantham_score_sim, grantham_score_combined := (grantham_score_mel)]
    # one species is NA and the other is SYN
    dt[!mel_pos & !sim_pos & (grantham_score_mel == 0 | grantham_score_sim == 0),
        grantham_score_combined := 0]

    # Summarize each group
    summarize_group <- function(classes, group_name) {
        group_dt <- dt[classification %in% classes]
        # combined mel/sim summary
        combined_summary <- group_dt[ grantham_score_combined > 0,
           .(   n_NS_snps = .N, 
                max_grantham = max(grantham_score_combined), 
                median_grantham = median(grantham_score_combined),
                mean_grantham = mean(grantham_score_combined)
           ), by = .(gene_id_fbgn_fromGFF, gene_name, longest_transcript_length)
        ][, `:=`(
                group = group_name,
                dataset = dataset_name,
                score_type = "Combined",
                snp_density = n_NS_snps / longest_transcript_length
            )
        ]
    # summary per species
    summarize_species <- function(score_col, species_name) {
            group_dt[
                get(score_col) > 0,
                .(
                    n_NS_snps = .N,
                    max_grantham = max(get(score_col)),
                    median_grantham = median(get(score_col)),
                    mean_grantham = mean(get(score_col))
                ),
                by = .(
                    gene_id_fbgn_fromGFF,
                    gene_name,
                    longest_transcript_length
                )
            ][
                ,
                `:=`(
                    group = group_name,
                    dataset = dataset_name,
                    score_type = species_name,
                    snp_density = n_NS_snps / longest_transcript_length
                )
            ]
        }

        mel_summary <- summarize_species("grantham_score_mel", "mel")
        sim_summary <- summarize_species("grantham_score_sim", "sim")

        rbindlist(list(combined_summary, mel_summary, sim_summary), use.names = TRUE, fill = TRUE)
    }

    rbindlist(
        list(
            summarize_group(NS_tsp, "TSP"),
            summarize_group(NS_conv, "Convergent"),
            summarize_group(NS_shared, "Shared")
        ), use.names = TRUE, fill = TRUE)
}

# All-SNP class definitions for the density numerator
all_snp_groups <- list(
    TSP = tsp_classes,          # A, B
    Convergent = conv_classes, # F, G, O, P, X, Y
    Shared = cand_classes       # A, B, F, G, O, P, X, Y
)

get_all_snp_density <- function(dt, dataset_name, gene_lengths) {
    dt <- copy(dt)
    dt <- merge(
        dt,
        gene_lengths[, .(gene_id_fbgn, longest_transcript_length)],
        by.x = "gene_id_fbgn_fromGFF",
        by.y = "gene_id_fbgn",
        all.x = TRUE
    )

    rbindlist(
        lapply(names(all_snp_groups), function(group_name) {
            classes <- all_snp_groups[[group_name]]
            dt[classification %in% classes &
                   !is.na(longest_transcript_length) &
                   longest_transcript_length > 0,
               .(
                   n_total_snps = .N,
                   total_snp_density = .N / longest_transcript_length
               ),
               by = .(gene_id_fbgn_fromGFF, gene_name, longest_transcript_length)
            ][, `:=`(
                    group = group_name,
                    dataset = dataset_name
                )
            ]
        }),
        use.names = TRUE,
        fill = TRUE
    )
}

sim_gff <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.gff"
mel_gff <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff"


poolseq_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/poolseq_dest/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.shared.bothMelSim.classed.rds")
inbred_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.11_2026.shared.bothMelSim.0.05MAFfilter.classed.rds")
gene_col <- "gene_id_fbgn_fromGFF"

cand_classes <- c("A", "B", "F", "G", "O", "P", "X", "Y")
tsp_classes <- c("A", "B")
conv_classes <- c("F", "G", "O", "P", "X", "Y")

NS_tsp <- c("A")
NS_conv <- c("F", "O", "X")
NS_shared <- c("A", "C", "E", "F", "H", "J", "K", "N", "O", "Q", "S", "T", "W", "X")

poolseq_dt <- add_grantham(poolseq_dt)
inbred_dt <- add_grantham(inbred_dt)


### max grantham score vs (total SNPs)/(length of longest transcript)
# like max grantham score of TSPs vs total TSPs/length transcript
# max grantham score of conv vs total conv / length transcript
# max grantham score of shared vs total shared / length transcript

# maybe also look at mean grantham score vs all those 

mel_gene_lengths <- get_longest_transcript_lengths(mel_gff)

poolseq_grantham <- summarize_grantham(
    poolseq_dt,
    dataset_name = "PoolSeq",
    gene_lengths = mel_gene_lengths
)

inbred_grantham <- summarize_grantham(
    inbred_dt,
    dataset_name = "Inbred",
    gene_lengths = mel_gene_lengths
)

grantham_gene_metrics <- rbind(
    poolseq_grantham,
    inbred_grantham,
    use.names = TRUE,
    fill = TRUE
)

# number excluded grantham scores 
poolseq_dt[
    grantham_score_mel > 0 &
    grantham_score_sim > 0 &
    grantham_score_mel != grantham_score_sim,
    .N
]

inbred_dt[
    grantham_score_mel > 0 &
    grantham_score_sim > 0 &
    grantham_score_mel != grantham_score_sim,
    .N
]

poolseq_dt[
    grantham_score_mel > 0 &
    grantham_score_sim > 0 &
    grantham_score_mel != grantham_score_sim,
    .N,
    by = classification
][order(-N)]
inbred_dt[
    grantham_score_mel > 0 &
    grantham_score_sim > 0 &
    grantham_score_mel != grantham_score_sim,
    .N,
    by = classification
][order(-N)]
############################################################
### PLOTS ###
combined_metrics <- grantham_gene_metrics[score_type == "Combined"]

# 1. max grantham score vs snp density 
plot_max_grantham <- ggplot(
    combined_metrics,
    aes(x = snp_density, y = max_grantham)
) +
    geom_point(alpha = 0.4, size = 1.2) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    scale_x_log10() +
    facet_grid(dataset ~ group) +
    theme_classic() +
    labs(
        x = "Nonsynonymous SNPs / longest CDS transcript length",
        y = "Maximum Grantham score",
        title = "Maximum amino acid substitution severity per gene"
    )

plot_max_grantham

# 2. mean grantham score vs snp density
plot_mean_grantham <- ggplot(
    combined_metrics,
    aes(x = snp_density, y = mean_grantham)
) +
    geom_point(alpha = 0.4, size = 1.2) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    scale_x_log10() +
    facet_grid(dataset ~ group) +
    theme_classic() +
    labs(
        x = "Nonsynonymous SNPs / longest CDS transcript length",
        y = "Mean Grantham score",
        title = "Mean amino acid substitution severity per gene"
    )

plot_mean_grantham

# 3. median grantham score vs snp density 
plot_median_grantham <- ggplot(
    combined_metrics,
    aes(x = snp_density, y = median_grantham)
) +
    geom_point(alpha = 0.4, size = 1.2) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    scale_x_log10() +
    facet_grid(dataset ~ group) +
    theme_classic() +
    labs(
        x = "Nonsynonymous SNPs / longest CDS transcript length",
        y = "Median Grantham score",
        title = "Median amino acid substitution severity per gene"
    )

plot_median_grantham



### quartiles of max grantham score vs snp density (number NS SNPs / length of longest CDS transcript) 
library(data.table)
library(ggplot2)

combined_metrics <- grantham_gene_metrics[score_type == "Combined"]
combined_metrics[, group := factor(    group,    levels = c("TSP", "Convergent", "Shared"))]
# Create density quintiles separately within each dataset
combined_metrics[, density_bin := factor(
    dplyr::ntile(snp_density, 4),
    levels = 1:4,
    labels = c("Q1", "Q2", "Q3", "Q4")),
    by = dataset]

make_grantham_bin_plot <- function(dt, dataset_title) {    
    ggplot(dt, aes(x = density_bin, 
    # y = max_grantham
    y = mean_grantham
    )) +        
    geom_boxplot(            
        width = 0.65,            
        outlier.shape = NA,            
        fill = "grey85",            
        color = "grey30"        
    ) +        
    geom_jitter(            
        width = 0.15,            
        alpha = 0.12,            
        size = 0.4        
    ) +        
    facet_wrap(~ group, nrow = 1) +        
    theme_classic() +        
    labs(            
        title = dataset_title,            
        x = "NS SNPs / longest CDS transcript length quartile (Q1 = lowest; Q4 = highest)",            
        # y = "Maximum Grantham score"        
        y = "Mean Grantham score"
    ) +        
    theme(            
        strip.background = element_blank(),            
        strip.text = element_text(size = 11),            
        axis.title = element_text(size = 10),            
        axis.text = element_text(size = 9),            
        plot.title = element_text(size = 12),            
        panel.spacing = unit(1, "lines")        
        )
}
plot_inbred_bins <- make_grantham_bin_plot(    
    combined_metrics[dataset == "Inbred"],    
    "Inbred")
plot_poolseq_bins <- make_grantham_bin_plot(    
        combined_metrics[dataset == "PoolSeq"],    
        "PoolSeq")

plot_inbred_bins
plot_poolseq_bins


### look at max grantham score vs (number all snps / length transcript) -- different from number NS snps / length transcript
### TEST: Mean Grantham vs ALL SNP density
### Grantham remains calculated on NS SNPs only

# Count ALL SNPs per gene, without filtering classification
poolseq_all_counts <- poolseq_dt[
    !is.na(gene_id_fbgn_fromGFF),
    .(n_all_snps = .N),
    by = gene_id_fbgn_fromGFF
]

inbred_all_counts <- inbred_dt[
    !is.na(gene_id_fbgn_fromGFF),
    .(n_all_snps = .N),
    by = gene_id_fbgn_fromGFF
]


# Start from existing NS-only Grantham metrics
test_metrics <- copy(
    grantham_gene_metrics[score_type == "Combined"]
)


# Merge all-SNP counts from the corresponding dataset
test_metrics[
    dataset == "PoolSeq",
    n_all_snps := poolseq_all_counts$n_all_snps[
        match(gene_id_fbgn_fromGFF,
              poolseq_all_counts$gene_id_fbgn_fromGFF)
    ]
]

test_metrics[
    dataset == "Inbred",
    n_all_snps := inbred_all_counts$n_all_snps[
        match(gene_id_fbgn_fromGFF,
              inbred_all_counts$gene_id_fbgn_fromGFF)
    ]
]


# Calculate ALL-SNP density using existing transcript length
test_metrics[
    ,
    total_snp_density := n_all_snps / longest_transcript_length
]


# Quartiles separately within each dataset (same approach as your original)
test_metrics[
    ,
    total_density_bin := factor(
        dplyr::ntile(total_snp_density, 4),
        levels = 1:4,
        labels = c("Q1", "Q2", "Q3", "Q4")
    ),
    by = dataset
]
test_metrics[, group := factor(
    group,
    levels = c("TSP", "Convergent", "Shared")
)]

# Plot function
make_total_density_plot <- function(dt, dataset_title) {

    ggplot(dt, aes(x = total_density_bin, y = mean_grantham)) +
        geom_boxplot(
            width = 0.65,
            outlier.shape = NA,
            fill = "grey85",
            color = "grey30"
        ) +
        geom_jitter(
            width = 0.15,
            alpha = 0.12,
            size = 0.4
        ) +
        facet_wrap(~group, nrow = 1) +
        theme_classic() +
        labs(
            title = dataset_title,
            x = "Number of SNPs / longest CDS transcript length quartile (Q1 = lowest; Q4 = highest)",
            y = "Mean Grantham score (NS SNPs only)"
        ) +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(size = 11),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 9),
            plot.title = element_text(size = 12),
            panel.spacing = unit(1, "lines")
        )
}



# Separate figures
plot_inbred_all_snp_bins <- make_total_density_plot(
    test_metrics[dataset == "Inbred"],
    "Inbred"
)

plot_poolseq_all_snp_bins <- make_total_density_plot(
    test_metrics[dataset == "PoolSeq"],
    "PoolSeq"
)

plot_inbred_all_snp_bins
plot_poolseq_all_snp_bins