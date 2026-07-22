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
out_csv <- "/scratch/ejy4bu/drosophila/GO/gowinda/gowindaRunsStats/gowindaStats_FDR0.05.csv"
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
file <- "/scratch/ejy4bu/drosophila/GO/gowinda/results/MAF1filter_polyAF/bg_speciesSpecific_noMAF/gowinda_XY_1_polyAF.txt"
results <- read.delim(file, header=FALSE, col.names=cols)
setDT(results)

results$Description <- go_file$name[results$GO.id]
results <- merge(results, go_dt, by="GO.id", all.x=T)

nrow(unique(results %>% filter(FDR < 0.05)))
nrow(unique(results %>% filter(p.value < 0.05)))

append_gowinda_summary(results, file, out_csv)



### FIGURES
library(data.table)
library(ggplot2)

csv <- fread(out_csv)

# class="AB"
# class="FGOPXY"
class="ABFGOPXY"
# class="XY"
# bg="speciesSpecific_noMAF"
bg="sharedOnly_noMAF"
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

### color code by ontology to easily examine BP GO terms
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

id <- "GO:0045879"

# data.table(
#     GO.id = id,
#     Name = Term(id),
#     Ontology = Ontology(id),
#     Definition = Definition(id)
# )


inspect_go <- function(go_ids){

 rbindlist(lapply(go_ids, function(go_id){
    genes <- unique(unlist(strsplit(gaf[GO.id == go_id, Gene.ids],"\\s+")))
    data.table(
        GO.id = go_id,
        Name = Term(go_id)[1],
        Ontology = Ontology(go_id)[1],
        N_genes = length(genes),
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