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

# results <- fread("/scratch/ejy4bu/drosophila/gowinda/results/gowinda_A_gene.txt", header=FALSE)

# results <- read.delim("/scratch/ejy4bu/drosophila/gowinda/results/gowinda_AB_gene.txt", header=FALSE, col.names=cols)
dir <- "/scratch/ejy4bu/drosophila/gowinda/maf_filter_mel5"
files_list <- list.files(path = dir, pattern = "\\.txt$", full.names = TRUE)

library(ontologyIndex)
go_file <- get_ontology("/scratch/ejy4bu/drosophila/gowinda/go.obo")


for (file_name in files_list) {

  group <- basename(file_name)
  group <- sub("^gowinda_", "", group)
  group <- sub("_snp_allBackground\\.txt$", "", group)

  assign(
    paste0("results_", group),
    read.delim(file_name, header = FALSE, col.names = cols)
  )
  file <- get(paste0("results_", group))
  file$Description <- go_file$name[file$GO.id]
  assign(paste0("results_", group), file)

  print(paste0(group, " with FDR<0.05: ", nrow(unique(file %>% filter(FDR < 0.05)))))
  rm(file)
}

#########################################################################
### GOdavid
#########################################################################

fwrite(
    data.table(gene_id = genes_AB),
    "/scratch/ejy4bu/drosophila/gowinda/GO_analysis/all_genes_AB.txt",
    col.names = FALSE
)

sig_AB <- unique(unlist(strsplit(
     results_AB[results_AB$FDR < 0.05, "Genes"],
     ","
 )))

#####################################################
# add ontology 
library(AnnotationDbi)
library(GO.db)
go_dt <- data.table(GO.id = keys(GO.db))
go_dt[, ontology := Ontology(GO.id)]
# gowinda_results from line 22, modify for each results file
gowinda_results <- merge(gowinda_results, go_dt, by="GO.id", all.x=T)

#####################################################
# slim
library(ontologyIndex)
library(data.table)
go_basic <- get_ontology("/scratch/ejy4bu/drosophila/gowinda/go-basic.obo", propagate_relationships = "is_a")
go_slim <- get_ontology("/scratch/ejy4bu/drosophila/gowinda/goslim_drosophila.obo")
slim_ids <- go_slim$id

# gowinda_results <- results_AB
gowinda_results <- results_FGOPXY
setDT(gowinda_results)

# get slim ancestors 
get_ancestors <- function(go_term, ontology, slim_id) {
  if (is.na(go_term) || !(go_term%in%names(ontology$ancestors))) return(NA_character_)
  
  ancestors <- ontology$ancestors[[go_term]]  # get all ancestors, also includes self term
  go_hits <- intersect(ancestors, slim_id)    # which ancestors are go slim terms
  if (length(go_hits) == 0) return(NA_character_)
  paste(go_hits, collapse=";")                # could be multiple terms
}

gowinda_results[, ancestors := lapply(GO.id, function(x) go_basic$ancestors[[x]])]
# gowinda_results[, slim_terms:= mapply(get_ancestors, GO.id, MoreArgs= list(ontology=go_basic, slim_id=slim_ids))]
slim_names <- data.table(
  slim_id = go_slim$id,
  slim_name = go_slim$name
)

# expand multi-ancestor terms
gowinda_slim <- gowinda_results[!is.na(slim_terms)] 
gowinda_slim <- gowinda_slim[, .(slim_id = unlist(strsplit(slim_terms, ";"))), by = setdiff(names(gowinda_slim), "slim_terms")]
gowinda_slim <- merge(gowinda_slim, slim_names, by="slim_id", all.x=T)

# # hypothesis GO terms ******* check these terms ******* (chosen from GO hierarchy)
# immune_root <- "GO:0006955"
# immune_roots <- c(
#     immune_root,
#     go_basic$children[[immune_root]]
#   )

# resistance_root <- "GO:0006805"  # xenobiotic metabolic process
# resistance_roots <- c(
#   resistance_root,
#   go_basic$children[[resistance_root]])


# # immune <- c(
# "GO:0002376", "GO:0045087", "GO:0045088",
# )
  # immune response
  # defense response
  # innate immune response


immune_roots <- slim_names[
  grepl("immune|defense|biotic", slim_name, ignore.case = TRUE)
  , slim_id
]
resistance_roots <- slim_names[
  grepl("xenobi|chemical|oxid|glutathione|drug", slim_name, ignore.case = TRUE)
  , slim_id
]

# resistance <- c()
  # response to chemical
  # ???

# immune_roots <- c(
#   "GO:0002376",
#   "GO:0006955",
#   "GO:0045087",
#   "GO:0009607"
# )

# resistance_roots <- c(
#   "GO:0009410",
#   "GO:0042221",
#   "GO:0055114",
#   "GO:0006629"
# )


# gowinda_slim[, hypothesis_category := fcase(
#   slim_id%in% immune, "immune",
#   slim_id%in% resistance, "resistance",
#   category=="biological_process", "other_BP",
#   category=="molecular_function", "other_MF",
#   default="other"
# )]

# gowinda_results[, ancestors := mapply(GO.id, get_ancestors, ontology = go_basic)]
# gowinda_results[, slim_list := strsplit(slim_terms, ";")]

# gowinda_results[, category := fcase(
#   sapply(slim_list, function(x) any(x %in% immune_roots)), "immune",
#   sapply(slim_list, function(x) any(x %in% resistance_roots)), "resistance",
#   default = "other"
# )]

gowinda_results[, category := fcase(
  sapply(ancestors, function(x) any(x %in% immune_roots)), "immune",
  sapply(ancestors, function(x) any(x %in% resistance_roots)), "resistance",
  default = "other"
)]

table(gowinda_results$category)

##############################################################
# trying functional classification
library(GO.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)

# gowinda_results <- results_AB
gowinda_results <- results_FGOPXY
go_ids <- unique(gowinda_results$GO.id)

go_map <- AnnotationDbi::select(
  GO.db,
  keys = go_ids,
  columns = c("TERM"),
  keytype = "GOID"
)

# go_map$category <- "Other"
# go_map$category[grepl("immune|defense|response to|virus|bacteria", go_map$TERM, ignore.case=TRUE)] <- "Immune/Defense"
go_map$category[grepl("metabolic|catabolic|biosynthetic|lipid|carbohydrate", go_map$TERM, ignore.case=TRUE)] <- "Metabolism"
go_map$category[grepl("transport|ion|transmembrane|localization", go_map$TERM, ignore.case=TRUE)] <- "Transport"
go_map$category[grepl("development|morphogenesis|pattern specification|differentiation", go_map$TERM, ignore.case=TRUE)] <- "Development"
go_map$category[grepl("signaling|signal transduction|receptor", go_map$TERM, ignore.case=TRUE)] <- "Signaling"
go_map$category[grepl("DNA|repair|replication|chromatin|cell cycle", go_map$TERM, ignore.case=TRUE)] <- "Genome maintenance"

# go_map$category <- "Other"
# IMMUNE
go_map$category[grepl("immune|defense|response to bacter|response to virus|antimicrobial", go_map$TERM, ignore.case = TRUE)] <- "Immune"
# ADAPTIVE / PESTICIDE / DETOX / RESISTANCE
go_map$category[grepl("detox|xenobiotic|drug|pesticide|insecticide|metabolic process|glutathione|cytochrome P450|oxidation|response to chemical",
  go_map$TERM,ignore.case = TRUE)] <- "Adaptive resistance"

plot <- as.data.frame(table(go_map$category))
plot2 <- as.data.frame(table(go_map$category[go_map$category%in%c("Immune", "Adaptive resistance")]))

ggplot(plot2, aes(x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat = "identity") +
  xlab("Functional Category") +
  ylab("Number of GO terms") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# junk

nrow(unique(results %>% filter(FDR < 0.05)))

### add GO Description to results file
library(ontologyIndex)
go_file <- get_ontology("/scratch/ejy4bu/drosophila/gowinda/go.obo")

head(go_file$name) # format: named vector GO:0000001 -> "mitochondrion inheritance"
results$Description <- go_file$name[results$GO.id]