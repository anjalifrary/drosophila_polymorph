library(data.table)
library(dplyr)


# column 1: the GO term
# column 2: on the average this number of genes are found per simulation for the given GO category. In --mode gene every gene is only counted once whereas in --mode snp a single gene may be counted several times dependent on the SNP
# column 3: using the candidate SNPs this number of genes was found for the given GO category. In --mode gene every gene is only counted once whereas in --mode snp a single gene may be counted several times dependent on the SNP
# column 4: p-value (uncorrected for multiple testing)
# column 5: FDR (p-value after adjustment for multiple testing)
# column 6: the number of genes (uniq) found for the given GO category
# column 7: the number of genes that could at most be found for the given GO category, i.e.: genes of the given GO category that have an corresponding entry in the annotation file and contain at least one SNP
# column 8: total number of genes for the given GO category in the GO association file
# column 9: description of the given GO term
# column 10: comma separated list of the gene_ids found for the given GO category

cols <- c("GO.id", "SimulatedGenes", "ObservedGenes", 
  "p.value", "FDR", 
  "N_GenesFound", "N_GenesMax", "N_GenesTotal",
  "Description", "Genes")



library(ontologyIndex)
go_file <- get_ontology("/scratch/ejy4bu/drosophila/GO/gowinda/go.obo")

#####################################################
# add ontology 
library(AnnotationDbi)
library(GO.db)
go_dt <- data.table(GO.id = keys(GO.db))
go_dt[, ontology := Ontology(GO.id)]
# gowinda_results from line 22, modify for each results file
# gowinda_results <- merge(gowinda_results, go_dt, by="GO.id", all.x=T)


parse_gowinda_metadata <- function(file){
    fname <- basename(file)
    dir_parts <- strsplit(dirname(file), "/")[[1]]

    # filename
    # gowinda_XY_1_polyAF.txt
    fname_parts <- strsplit(sub("\\.txt$", "", fname), "_")[[1]]

    classes <- fname_parts[2]
    maf_value <- fname_parts[3]
    maf_def <- fname_parts[4]

    # directory pieces
    maf_dir <- dir_parts[grepl("MAF", dir_parts)]
    bg_dir <- dir_parts[grepl("^bg_", dir_parts)]

    background <- sub("^bg_", "", bg_dir)

    data.table(
        classes = classes,
        MAF_value = as.numeric(maf_value),
        MAF_def = maf_def,
        background = background
    )
}

append_gowinda_summary <- function(results, file, csv_file, filter_col="FDR", threshold=0.05){
    meta <- parse_gowinda_metadata(file)
    sig <- unique(results[get(filter_col) < threshold, GO.id])
    new_row <- cbind(
        meta,
        data.table(
            threshold = paste0(filter_col,"<",threshold),
            N_GOTerms = length(sig),
            GO.ids = if(length(sig)) paste(sig, collapse=";") else NA_character_
        )
    )

    if(file.exists(csv_file)){
        out <- fread(csv_file)

        # check if this exact analysis already exists
        duplicate <- out[
            classes == new_row$classes &
            MAF_value == new_row$MAF_value &
            MAF_def == new_row$MAF_def &
            background == new_row$background &
            threshold == new_row$threshold
        ]

        if(nrow(duplicate) > 0){
            message("Already exists: skipping ", basename(file))
            return(invisible(out))
        }
        out <- rbind(out, new_row, fill=TRUE)

    } else {
        out <- new_row
        file.create(csv_file)
    }

    fwrite(out, csv_file)
}


# csv headers
# MAF_value | MAF_def | background | classes | statThreshold | N_GOTerms | GO.ids
out_csv <- "/scratch/ejy4bu/drosophila/GO/gowinda/gowindaRunsStats/gowindaStats.8-6-2026.csv"
dir.create("/scratch/ejy4bu/drosophila/GO/gowinda/gowindaRunsStats/")

### to loop over a single results directory:

dir <- "/scratch/ejy4bu/drosophila/GO/gowinda/results/"
files_list <- list.files(path = dir, pattern="gowinda_.*txt", recursive = TRUE, full.names = TRUE)

for (file_name in files_list) {
    message("Processing: ", file_name)

    # skip empty files
    if(file.info(file_name)$size == 0){
        message("Skipping empty file: ", file_name)
        next
    }
    results <- fread(file_name, header=FALSE, col.names=cols)
    append_gowinda_summary(results, file_name, out_csv)
}


### to inspect a single file:
file <- "/scratch/ejy4bu/drosophila/GO/gowinda/results/MAF1filter_globalAF/bg_speciesSpecific_MAF/gowinda_AB_1_globalAF.txt"
results <- read.delim(file, header=FALSE, col.names=cols)
setDT(results)

results$Description <- go_file$name[results$GO.id]
results <- merge(results, go_dt, by="GO.id", all.x=T)

nrow(unique(results %>% filter(FDR < 0.05)))
nrow(unique(results %>% filter(p.value < 0.05)))

append_gowinda_summary(results, file, out_csv)

################################################################################

### FIGURES
library(data.table)
library(ggplot2)

csv <- fread(out_csv)

class="AB"
# class="FGOPXY"
# class="ABFGOPXY"
# class="XY"
# bg="speciesSpecific_MAF"
# bg="sharedOnly_MAF"
bg="melOnly_MAF"
maf_def="polyAF"
# maf_def="globalAF"

csv <- as.data.table(csv)

plot_counts <- csv[ 
    classes==class & 
    background==bg & 
    MAF_def==maf_def,
    .(classes, MAF_value, MAF_def, background, threshold, N_GOTerms, GO.ids)]


### Persistance of specific GO terms:
go_long <- plot_counts[N_GOTerms > 0, 
    .(GO.id = unlist(strsplit(GO.ids, ";"))),
    by = .(MAF_value)
]
go_long[, persistence := uniqueN(MAF_value), by = GO.id]

# y = number significant go terms
# x = MAF 
ggplot(plot_counts, aes(x=MAF_value, y=N_GOTerms)) +
    geom_point(size=3) +
    geom_line() +
    scale_x_continuous(breaks=unique(plot_counts$MAF_value)) +
    labs(
        x="MAF filter (%)",
        y="Number of significant GO terms", 
        title=paste0(class, " , ", bg, " background")
    ) +
    theme_classic()

### add color code by ontology to easily examine BP GO terms
ggplot(go_long,
       aes(x = factor(MAF_value),
           y =  reorder(GO.id, -persistence))) +
    geom_tile() +
    labs(
        x = "MAF filter (%)",
        y = "GO term",
        title=paste0(class, " , ", bg, " background")

    ) +
    theme_classic() + 
    theme(axis.text.y = element_text(size = 5)) + 
    scale_y_discrete(limits = rev)


### inspect go ids
gaf <- fread(
    "/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/flybase_gaf_go.txt",
    header=FALSE, col.names=c("GO.id", "GO.id", "Gene.ids")
)

id <- c("GO:0005615", "GO:0005886", "GO:0032590")

# data.table(
#     GO.id = id,
#     Name = Term(id),
#     Ontology = Ontology(id), 
#     Definition = Definition(id)
# )

# returns table of 1 row per go id 

inspect_go <- function(go_ids){
 rbindlist(lapply(go_ids, function(go_id){
    genes <- unique(unlist(strsplit(gaf[GO.id == go_id, Gene.ids],"\\s+")))
    data.table(
        GO.id = go_id,
        Name = Term(go_id)[1],
        Ontology = Ontology(go_id)[1],
        N_genes_possible = length(genes),
        Persistence = unique(go_long[GO.id==go_id, persistence])[1],
        Definition = Definition(go_id)[1],
        Genes = list(genes)
    )
}))}

View(inspect_go(id))
### one row per gene
inspect_go_genePerRow <- function(go_id){

    genes <- unique(unlist(strsplit(gaf[GO.id == go_id, Gene.ids],"\\s+")))

    data.table(
        GO.id = go_id,
        Name = Term(go_id),
        Ontology = Ontology(go_id),
        Definition = Definition(go_id),
        Gene.id = genes
    )
}

View(inspect_go_genePerRow(id))


### extract from table used to make figures
sort(unique(go_long$GO.id))

# sort by persistence
unique(go_long[, .(GO.id, persistence)])[order(-persistence, GO.id)]

persistent <- unique(
    go_long[persistence >= 3, GO.id]
)

persistent_info <- inspect_go(persistent)
View(persistent_info)

all_info <- inspect_go(go_long$GO.id)
View(all_info)
#####################################################################################

### 7/29 build long table to keep genes per GO term per threshold

# keep parse_gowinda_metadata() from above

build_long_gowinda <- function(files, filter_col = "FDR", threshold= 0.05) {
    rbindlist(lapply(files, function(file) {
        if (file.info(file)$size == 0 ){
            message("skipping empty file: ", file)
            return(NULL)
        }
        results <- fread(file, header=FALSE, col.names = cols)
        meta <- parse_gowinda_metadata(file)
        results[, Description := go_file$name[GO.id]]
        results <- merge(results,go_dt,by = "GO.id",all.x = TRUE)
        results[, Definition := Definition(GO.id)]
        results[, GeneListFound := strsplit(Genes, ",")]
        meta[, file := basename(file)]
        sig <- results[get(filter_col) < threshold]
        if(nrow(sig) == 0) return(NULL)
        cbind(meta, sig[, .(GO.id, Description, ontology, Definition, FDR, p.value, N_GenesFound, N_GenesMax, N_GenesTotal, GeneListFound)])
    }), fill = T)
}

maf_inputs <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.49)

# after building long table: same plots as before
plot_go_counts <- function(dt, class_name, bg_name, maf_def_name) {
    counts <- dt[classes == class_name & background == bg_name & MAF_def == maf_def_name,
                        .(N_GOTerms = uniqueN(GO.id)), , by = MAF_value]
    counts <- merge(data.table(MAF_value = maf_inputs * 100), counts, by = "MAF_value", all.x = TRUE)
    counts[is.na(N_GOTerms), N_GOTerms := 0]

    ggplot(counts, aes(x = MAF_value, y = N_GOTerms)) +
    geom_point(size = 3) + geom_line() +
    scale_x_continuous(limits = c(0, 50), breaks = maf_inputs*100) +
        # breaks = seq(0, 50, 10))+
        # breaks = unique(counts$MAF_value)) +
    labs(x = "MAF filter (%)", y = "Number of significant GO terms",
         title = paste0(class_name, ", ", bg_name)) +
    theme_classic()
}

plot_go_persistence_heatmap <- function(dt, class_name, bg_name, maf_def_name) {
  sub <- dt[classes == class_name & background == bg_name & MAF_def == maf_def_name]
  sub[, persistence := uniqueN(MAF_value), by = GO.id]
  ggplot(sub, aes(x = factor(MAF_value), y = reorder(GO.id, -persistence))) +
    geom_tile() +
    labs(x = "MAF filter (%)", y = "GO term",
         title = paste0(class_name, ", ", bg_name)) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 5)) +
    scale_y_discrete(limits = rev)
}

# # genes enriched for 1 GO term at each MAF threshold to see if same genes drive signal or if they vary across MAF threshold
# genes_by_threshold <- function(dt, go_id) {
#   rows <- copy(dt[GO.id == go_id, .(classes, background, MAF_def, MAF_value, GeneList)])
#   if(nrow(rows) == 0) {
#     message("GO not found: ", go_id)
#     return(null)
#   }

#   setorder(rows, MAF_value)
# }

# view which genes persist across thresholds:
gene_persistence_matrix <- function(dt, go_id, class_name, bg_name, maf_def_name) {
  rows <- copy(dt[GO.id == go_id & classes == class_name & background == bg_name & MAF_def == maf_def_name])
  if(nrow(rows) == 0) {
    message("GO not found: ", go_id)
    return(null)
  }

  setorder(rows, MAF_value)

  gene_sets <- rows$GeneList
  names(gene_sets) <- rows$MAF_value
 
  all_genes <- sort(unique(unlist(gene_sets)))
  presence <- sapply(gene_sets, function(g) all_genes %in% g)
  out <- data.table(Gene.id = all_genes,presence)
  setnames(out, old = names(gene_sets), new = paste0("MAF_" , names(gene_sets)))
  
  out
}

# to use:
dir <- "/scratch/ejy4bu/drosophila/GO/gowinda/results/"
files_list <- list.files(path = dir, pattern="gowinda_.*txt", recursive = TRUE, full.names = TRUE)

long_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/gowinda_results_all_longFormat.rds")
long_dt <- build_long_gowinda(files_list)
# saveRDS(long_dt, "/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/gowinda_results_all_longFormat.rds")


View(long_dt)
View(long_dt[classes=="AB" & background=="speciesSpecific_MAF" & MAF_def == "polyAF", ])

plot_go_counts(long_dt, "AB", "speciesSpecific_MAF", "polyAF")
plot_go_persistence_heatmap(long_dt, "AB", "speciesSpecific_MAF", "polyAF")

View(genes_by_threshold(long_dt, "GO:0006122"))
View(gene_persistence_matrix(long_dt, "GO:0006122", "AB", "speciesSpecific_MAF", "polyAF"))
#####################################################################################

### 7/29 attempt at consolidating GO analysis / visualization (see GO_SemanticSimilarity.R and gowinda_analysis.R for more)
library(data.table)
library(ggplot2)
library(GOSemSim)
library(AnnotationDbi)
library(GO.db)
library(org.Dm.eg.db)

## semantic similarity

# set up:
semdata <- function(ont = "BP", annoDb = "org.Dm.eg.db") {
    godata(annoDb = annoDb, ont = ont)
}

# function to calculate pairwise similarity + hierarch. clustering + gene ann
# for 1 set of GO ids
go_cluster_summary <- function(go_ids, semData, 
    gaf_file="/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/flybase_gaf_go.txt",
    k = 10, measure = "Wang") {
  
    go_ids <- unique(go_ids)
    go_ids <- go_ids[!is.na(go_ids)]

    if(length(go_ids) < 2){
        message("Need at least 2 GO terms")
        return(NULL)
    }
                
    similarity <- mgoSim(go_ids, go_ids, semData=semData, measure=measure, combine=NULL)

    hc <- hclust(as.dist(1-similarity))

    clusters <- cutree(hc, k=k)

    gaf <- fread(gaf_file, header = FALSE, col.names = c("GO.id", "GO.id2", "Gene.ids"))

    info <- rbindlist(lapply(go_ids, function(id) {                                   
        genes <- unique(unlist(strsplit(gaf[GO.id == id, Gene.ids], "\\s+")))
        data.table(
            GO.id = id, 
            name = Term(id)[1], 
            Ontology = Ontology(id)[1],
            N_genes_possible = length(genes),
            Definition = Definition(id)[1],
            GeneList = list(genes),
            N_genes_possible = length(genes)
        )
    }))   
    info[, cluster := clusters[GO.id]]
    setorder(info, cluster)
    list(sim_matrix = similarity, hclust = hc, info = info)

}

# function to compare similarity between 2 diff GO sets (like AB vs FGOPXY; AB maf1 vs maf2)
# combine BMA = collapse to 1 score; combine = NULL give full matrix on a per-term level
compare_diff_go_sets <- function(go_set1, go_set2, semData = semData, measure = "Wang", combine = "BMA") {
    mgoSim(unique(go_set1), unique(go_set2), semData = semData, measure = measure, combine = combine)
}


# example to run functions defined above:
dmGO <- semdata("BP")

dir <- "/scratch/ejy4bu/drosophila/GO/gowinda/results/"
files_list <- list.files(path = dir, pattern="gowinda_.*txt", recursive = TRUE, full.names = TRUE)

long_dt <- build_long_gowinda(files_list)
AB_terms <- unique(long_dt[classes=="AB" & background=="speciesSpecific_MAF" & MAF_def == "polyAF", GO.id])

res <- go_cluster_summary(AB_terms, dmGO)
View(res$info[order(cluster)])
plot(res$hc, cex=0.4)

AB_terms2 <- unique(long_dt[classes=="AB" & background=="sharedOnly_MAF" & MAF_def == "polyAF", GO.id])
compare_diff_go_sets(AB_terms, AB_terms2, dmGO)


########## in progress: 

library(data.table)

go <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/gowinda_results_all_longFormat.rds")
dat <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/noMAFfilter/subset_qualVar_ofInterest_7-20-2026.rds")

load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_simulans.17_06_2026.nlpTable.Rdata")
sim_nlp <- nlp 
rm(nlp)
sim_nlp[, poly_af_sim := poly_af]
sim_nlp[, global_af_sim := global_af]

load("/project/berglandlab/multispecies_endemism/data/collectiveAnalysis_version3/Drosophila_melanogaster.10_03_2026.nlpTable.paramask.genmap.busco.Rdata")
mel_nlp <- nlp
rm(nlp)   
mel_nlp[, poly_af_mel := poly_af]
mel_nlp[, global_af_mel := global_af]

mel_af <- mel_nlp[, .(chr, pos, poly_af_mel, global_af_mel)]
dat <- merge(dat, mel_af, by=c("chr", "pos"), all.x=T)

sim_af <- sim_nlp[, .(chr, pos, poly_af_sim, global_af_sim)]
dat <- merge(dat, sim_af, by=c("chr", "pos"), all.x=T)

candidates <- eval(go[classes=="AB"][MAF_value==10][MAF_def=="polyAF"][background=="speciesSpecific_noMAF"][order(FDR)][23,]$GeneListFound)[[1]]
candidates <- gsub("fbgn", "FBgn", candidates)

View(dat[gene_id_mel%in%candidates][effect_mel%like%"syn"][poly_af_mel>.1])
View(dat[gene_id_mel%in%candidates][effect_mel%like%"missense"][poly_af_mel>.1])

View(dat[gene_id_mel%in%candidates][effect_mel%like%"syn"][poly_af_sim>.1])
View(dat[gene_id_mel%in%candidates][effect_mel%like%"missense"][poly_af_sim>.1])

###########################################################################

### 8/6/2026:
library(ggplot2)
library(data.table)
library(ggnewscale)


long_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/gowinda_results_all_longFormat.rds")

### 1. inspect gene list and GO terms across different backgrounds / MAF def but within same group
# let's start with polyAF 
# melOnly_MAF; speciesSpecfic_MAF; sharedOnly_MAF

# group <- "AB"
group <- "ABFGOPXY"
# group <- "FGOPXY"
maf_def <- "polyAF"
bg_list <- c("melOnly_MAF", "speciesSpecific_MAF", "sharedOnly_MAF")

go_subset <- long_dt[MAF_def==maf_def]

maf_values <- c(0.5, 1, 2, 5, 10, 15, 20, 25, 30, 40, 49)

go_terms <- go_subset[
    classes==group & MAF_value%in%maf_values & background%in%bg_list, .(GO.id, Description, ontology, MAF_value, background, GeneListFound)
]

go_terms

go_terms_perBG <- dcast(
    unique(go_terms[,.(GO.id, background)]),
    GO.id ~ background,
    fun.aggregate = length,
    value.var = "background"
)

go_plot <- melt(go_terms_perBG, id.vars = "GO.id")
go_order <- go_plot[, .(persistence = sum(value)), by = GO.id][order(-persistence), GO.id]
go_plot[, GO.id := factor(GO.id, levels = go_order)]

ontology_plot <- unique(go_terms[, .(GO.id, ontology)])
ontology_plot[, GO.id := factor(GO.id, levels = go_order)]

ggplot(go_plot, aes(x=variable, y=GO.id)) + 
    geom_tile(aes(fill=value)) + 
    scale_fill_gradient(low = "white", high = "black", guide="none") + 
    ggnewscale::new_scale_fill() +
    geom_tile(data = ontology_plot,aes(x="Ontology", y=GO.id, fill=ontology)) +
    scale_fill_manual(values=c(
        "BP"="darkgreen",
        "MF"="steelblue",
        "CC"="darkred"), name="Ontology") +
    # geom_text(data=ontology_plot, aes(x=6, label=ontology), size=1.5) + 
    theme_classic() +
    labs(x="Background", y = "GO term") + 
    theme(
        axis.text.x = element_text(size = 6, angle=20, hjust=1),
        axis.text.y = element_text(size = 6), 
        axis.title = element_text(size = 6.5)
    ) 

### repeat for genes 
### and calculate semantic similarity metric 

### 2. inspect gene list and GO terms within same background / MAF def but across diff group 

bg = "speciesSpecific_MAF"
# bg = "melOnly_MAF"
# bg = "sharedOnly_MAF"
groups = c("AB", "FGOPXY", "ABFGOPXY")
maf_def = "polyAF"

### 3. merge on MKish stats to look at specific genes that are significant in both terms



### 4. looking at defense GO: GO.id==0006952
long_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/gowinda_results_all_longFormat.rds")

defense <- long_dt[GO.id=="GO:0006952", ]
defense_genes <- as.data.table(unique(defense[, unlist(GeneListFound)]))
setnames(defense_genes, "gene_id")

defense_genes[, gene_id := sub("^fbgn", "FBgn", gene_id)]

library(AnnotationDbi)
library(org.Dm.eg.db)

gene_info <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys = defense_genes$gene_id,
    keytype = "FLYBASE",
    columns = c("FLYBASE", "FLYBASECG", "SYMBOL", "GENENAME", "ENTREZID")
)

gene_info <- setnames(gene_info, "FLYBASE", "gene_id")

### add BUSCO data 
load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_melanogaster.11_08_2026.nlpTable.paramask.genmap.busco.repeatmask.wmdust.Rdata")
mel_nlp <- as.data.table(nlp)
rm(nlp)

# collapse BUSCO to 1 row per gene
mel_busco <- mel_nlp[ , .( BUSCO = if ( all(is.na(busco)) ) { NA_character_ } else { unique(na.omit(busco))[1] } ), by = gene ]

gene_info <- merge( gene_info, mel_busco, by.x = "FLYBASECG", by.y = "gene", all.x = TRUE )
defense_genes <- merge(defense_genes, gene_info, by="gene_id")

### add MKish stats 
MKish <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/adaptedMK/new_asymptotic_MK_longResults_polyAF_speciesSpecificBG.rds")

def_genes <- merge(defense_genes, MKish, by.x="gene_id", by.y="gene", all.x=T)

### merge on master candidate snps 
# masterCandidates <- readRDS("/scratch/ejy4bu/drosophila/GO/gowinda/candidateFiles/masterCandidateFile.rds")
# all_qualVar <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/noMAFfilter/all_quality_variants_clean.rds") # not classified
voi <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/noMAFfilter/voi_qualVar_ofInterest_07-20-2026.rds")
subset <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/noMAFfilter/subset_qualVar_ofInterest_final.rds") # this one has more matching snps
subset2 <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed.rds")
snps <- subset2
def_genelist <- unique(def_genes[, gene_id])
defense_snps <- snps[gene_id_mel%in%def_genelist, ]

defense_snps <- merge(defense_snps, mel_nlp[, .(chr, pos, global_af, poly_af, nLocales_poly, busco)], by=c("chr", "pos"))
defense_snps <- merge(defense_snps, defense_genes[, .(gene_id, SYMBOL, GENENAME)], by.x="gene_id_mel", by.y="gene_id")

# reorder cols
library(dplyr)
defense_snps <- defense_snps %>% relocate(gene_id_mel, .after=classification)
defense_snps <- defense_snps %>% relocate(SYMBOL, .after=gene_id_mel)
defense_snps <- defense_snps %>% relocate(GENENAME, .after=gene_id_mel)
defense_snps <- defense_snps %>% relocate(GENENAME, .after=SYMBOL)
defense_snps <- defense_snps %>% relocate(nLocales_poly, .after=aa_alt_sim)
defense_snps <- defense_snps %>% relocate(global_af, .after=nLocales_poly)
defense_snps <- defense_snps %>% relocate(poly_af, .after=nLocales_poly)

View(defense_snps[SYMBOL%like%"Drsl", ])

### i'm rly confused and there is a snp in Drsl3 that is in the same codon as another snp so what happened there??? 

### get snp counts in mel_only, sim_only, shared
### add AB counts, FGOPXY counts

### other data 
# number of candidate SNPs
# number of shared polymorphisms
# number of species-specific polymorphisms
# number of SNPs at each MAF threshold
# mean/median MAF
# maximum MAF
# number of synonymous vs nonsynonymous variants
# number of variants per codon
# potentially the number of variants classified as TSP candidates

### 5. looking at busco genes

load("/scratch/ejy4bu/drosophila/GO/gowinda/Drosophila_melanogaster.11_08_2026.nlpTable.paramask.genmap.busco.repeatmask.wmdust.Rdata")

long_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/gowinda_results_all_longFormat.rds")

enriched_genes <- unique(unlist(long_dt$GeneListFound))

enriched_genes <- sub("^fbgn", "FBgn", enriched_genes)

head(enriched_genes)

nlp_busco <- nlp[, .(
    busco = if (all(is.na(busco))) NA_character_ else unique(na.omit(busco))[1]
), by = gene]


gene_map <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys = enriched_genes,
    keytype = "FLYBASE",
    columns = c("FLYBASE", "FLYBASECG")
)


gene_map_busco <- merge(
    gene_map,
    nlp_busco,
    by.x = "FLYBASECG",
    by.y = "gene",
    all.x = TRUE
)
busco_genes <- unique(
    gene_map_busco[["FLYBASE"]][!is.na(gene_map_busco[["busco"]])]
)
long_dt[, n_BUSCOgenes := vapply(
    GeneListFound,
    function(x) {
        genes <- unique(sub("^fbgn", "FBgn", unlist(x)))
        sum(genes %in% busco_genes)
    },
    integer(1)
)]

saveRDS(long_dt, "/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/gowinda_results_all_longFormat.busco.rds")

# repeat for prop of busco genes (n_BUSCO / n_genes found)
library(stringr)

plot_dt <- long_dt[
    classes == "AB" &
    MAF_def == "polyAF" &
    background == "sharedOnly_MAF" &
    FDR < 0.05
]

plot_dt[, prop_BUSCO := n_BUSCOgenes / N_GenesFound]

plot_dt <- plot_dt[order(prop_BUSCO)]

ggplot(
    plot_dt,
    aes(
        x = prop_BUSCO,
        y = reorder(Description, prop_BUSCO)
    )
) +
    geom_point(size = 2.5) +
    scale_y_discrete(
        labels = \(x) stringr::str_wrap(x, width = 45)
    ) +
    theme_classic() +
    labs(
        x = "Prop of BUSCO genes",
        y = "GO term"
    ) +
    theme(
        axis.text.y = element_text(size = 4)
    )