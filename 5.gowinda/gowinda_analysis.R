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
dir <- "/scratch/ejy4bu/drosophila/gowinda/results/final"
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

  print(paste0(group, ": ", nrow(unique(file %>% filter(FDR < 0.05)))))
  rm(file)
}



nrow(unique(results %>% filter(FDR < 0.05)))

### add GO Description to results file
library(ontologyIndex)
go_file <- get_ontology("/scratch/ejy4bu/drosophila/gowinda/go.obo")

head(go_file$name) # format: named vector GO:0000001 -> "mitochondrion inheritance"
results$Description <- go_file$name[results$GO.id]