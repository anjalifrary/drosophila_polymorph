library(data.table)
library(ggplot2)


out_dir <- "/scratch/ejy4bu/drosophila/gds_analysis/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

in_rds <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variantsOfInterest.rds"
dt <- readRDS(in_rds)

dt[, species := fifelse(
    !is.na(ref_sim), "sim",
    fifelse(!is.na(ref_mel), "mel", NA_character_)
)]

dt <- dt[!is.na(species)]
dt <- dt[order(chr, pos)]

chroms <- unique(dt$chr)
# chroms <- chroms[1]

for (chrom in chroms) { 
    dt_chr <- dt[chr == chrom] # subset dt for chromosome
    chr_len <- max(dt[chr == chrom]$pos) # get chr length 
    wrap_size <- (chr_len / 4)+1 # divide the chromosome into 4 lines for wrapping
    dt_chr[, line := floor(pos / wrap_size)] # set length of each line 

    ##### Scatter plot #####
    p <- ggplot(dt_chr, aes(x=pos, y=species, color=species)) +
        geom_point(alpha = 0.4, size = 0.6, position = position_jitterdodge(jitter.height = 0.3)) +
        theme_bw() + 
        labs(
            title=paste("Chromosome ", chrom),
            x = "Genomic Position", 
            color = "Species"
        ) + 
        facet_wrap(~line, scales = "free_x", ncol = 1) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 6))

    ggsave(
    filename = file.path(out_dir, paste0("variant_wrapped_scatter_per_chr_", chrom, ".png")),
    plot = p,
    width = 12,
    height = 6,
    dpi = 300
)
} 

