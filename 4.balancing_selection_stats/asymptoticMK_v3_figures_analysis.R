# library(SeqArray)
library(data.table)
# library(doMC)
# registerDoMC(16)
library(ggplot2)

results <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/adaptedMK/new_asymptotic_MK_longResults_polyAF_speciesSpecificBG.rds")

nrow(results[!is.na(alpha)])

summary(results[, alpha])
summary(results[, logOR])

### FIGURES ###

### 1. histogram, breaks=100 
# hist symmetry should shift to center around 0 as MAF increases... ? 
# add -inf / inf bar on the ends ? 
maf = 0.005
group = "TSP"
MAF_def = "polyAF"

maf_inputs <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.49)
maf_inputs <- c(0.005, 0.01)
# maf_inputs <- c(0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.49)
for (maf in maf_inputs) {    
    maf_label = maf * 100
    subset <- results[MAF==maf & !is.na(logOR), ]
    par(mfrow = c(1,3))
    xlim <- range(subset$logOR)
    # ylim <- max(table(cut(subset$logOR, breaks=100)))

    for(g in c("TSP","CONV","BOTH")){
        vals <- subset[group == g, logOR]

        if(length(vals) > 0){
            hist(
                subset[group == g, logOR],
                breaks = 100,
                xlim = xlim,
                main = paste0(g, "_", maf_label, "% MAF"),
                xlab = "log(OR)"
            )   
        } else {
            plot.new()
            title(main = paste0(g, "_", maf_label, "% MAF\n(no data)"))
        }
    }
}
    # hist(subset[group=="TSP", alpha], breaks=100, main=paste0("TSP_", maf_label))
    # hist(subset[group=="CONV", alpha], breaks=100)
    # hist(subset[group=="BOTH", alpha], breaks=100)



### 2. box-whisker plots, alpha vs MAF with annotated label# genes for which alpha !is.na
counts <-results[!is.na(alpha), .N, by=.(MAF, group)]

ggplot(results[!is.na(alpha), ], aes(x=factor(MAF), y=alpha)) +
    geom_boxplot() +
    facet_wrap(~group) +
    theme_classic() +
    labs(x="MAF threshold", y="alpha") + 
    geom_text(data = counts,
        aes(x=factor(MAF), y=Inf, label=N),
        vjust=1.5,
        size=2
    ) + 
    theme(axis.text.x = element_text(size=7, angle=45, hjust=1))

### 2b. box-whisker log(OR) vs MAF
counts <-results[!is.na(logOR), .N, by=.(MAF, group)]

ggplot(results[!is.na(logOR), ], aes(x=factor(MAF), y=logOR)) +
    geom_boxplot() +
    facet_wrap(~group) +
    theme_classic() +
    labs(x="MAF threshold", y="log(OR)") + 
    geom_text(data = counts,
        aes(x=factor(MAF), y=Inf, label=N),
        vjust=1.5,
        size=2
    ) + 
    theme(axis.text.x = element_text(size=7, angle=45, hjust=1))



### 3. line graph alpha vs MAF where each line is a gene (like asymptotic application)

ggplot(results[!is.na(alpha)], 
       aes(x=MAF, y=alpha, group=gene)) +
    geom_line(alpha=0.1) +
    facet_wrap(~group) +
    theme_classic() +
    labs(x="MAF threshold", y="alpha")


gene_counts <- results[!is.na(alpha), .N, by=gene]
filter <- c(0, 5, 10, 15, 20)
for (threshold in filter) {
    stable_genes <- gene_counts[N >= threshold, gene]
    p <- ggplot(results[gene %in% stable_genes & !is.na(alpha)],
        aes(x=MAF, y=alpha, group=gene)) +
        geom_line(alpha=0.1) +
        facet_wrap(~group) +
        theme_classic() +
        labs(x="MAF threshold", y="alpha", title=paste0("Genes in >", threshold, " rows"))
    
    print(p)
}

### 4. line graph log(OR) vs MAF where each line is a gene (like asymptotic application)
ggplot(results[!is.na(logOR)], 
       aes(x=MAF, y=logOR, group=gene)) +
    geom_line(alpha=0.1) +
    facet_wrap(~group) +
    theme_classic() +
    labs(x="MAF threshold", y="log(OR)")


gene_counts <- results[!is.na(logOR), .N, by=gene]
filter <- c(0, 5, 10, 15, 20)
for (threshold in filter) {
    stable_genes <- gene_counts[N >= threshold, gene]
    p <- ggplot(results[gene %in% stable_genes & !is.na(logOR)],
        aes(x=MAF, y=logOR, group=gene)) +
        geom_line(alpha=0.1) +
        facet_wrap(~group) +
        theme_classic() +
        labs(x="MAF threshold", y="log(OR)", title=paste0("Genes in >", threshold, " rows"))
    
    print(p)
}

### 5. median trend lines 

ggplot(results[!is.na(alpha)], aes(x=MAF, y=alpha, group=group)) +
    geom_line(
        data=results[!is.na(alpha),
             .(alpha=median(alpha)),
             by=.(MAF,group)]
    ) +
    facet_wrap(~group) +
    theme_classic() +
    labs(x="MAF", y="median alpha")

ggplot(results[!is.na(logOR)], aes(x=MAF, y=logOR, group=group)) +
    geom_line(
        data=results[!is.na(logOR),
             .(logOR=median(logOR)),
             by=.(MAF,group)]
    ) +
    facet_wrap(~group) +
    theme_classic() + 
    labs(x="MAF", y="median log(OR)")