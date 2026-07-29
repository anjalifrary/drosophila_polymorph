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
out_csv <- "/scratch/ejy4bu/drosophila/GO/gowinda/gowindaRunsStats/gowindaStats.csv"
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
bg="speciesSpecific_MAF"
# bg="sharedOnly_MAF"
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

# after building long table: same plots as before
plot_go_counts <- function(dt, class_name, bg_name, maf_def_name) {
  counts <- dt[classes == class_name & background == bg_name & MAF_def == maf_def_name,
                      .N, by = MAF_value]
    ggplot(counts, aes(x = MAF_value, y = N)) +
    geom_point(size = 3) + geom_line() +
    scale_x_continuous(breaks = unique(counts$MAF_value)) +
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




# analyze_gowinda <- function(results, fdr = 0.05, semData = dmGO, ontology = "BP") {
#     results <- as.data.table(results)

#     sig <- results[
#         ontology == ontology &
#         FDR < fdr
#     ]

#     if (nrow(sig) == 0)
#         return(NULL)

#     ## parse gene list
#     sig[, GeneList := strsplit(Genes, ",")]

#     ## semantic similarity
#     sim <- mgoSim(
#         sig$GO.id,
#         sig$GO.id,
#         semData = semData,
#         measure = "Wang",
#         combine = NULL
#     )

#     hc <- hclust(as.dist(1 - sim))

#     sig[, cluster := cutree(hc, k = 10)]

#     list(
#         results = sig,
#         similarity = sim,
#         tree = hc
#     )
# }

# AB <- analyze_gowinda(results_AB)

# View(AB$results)

# heatmap(AB$similarity)

# plot(AB$tree)