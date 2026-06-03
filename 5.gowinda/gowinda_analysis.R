library(data.table)
library(dplyr)


cols <- c("GO.id", "SimulatedGenes", "ObservedGenes", 
  "p.value", "FDR", 
  "N_GenesFound", "N_GenesMax", "N_GenesTotal",
  "Description", "Genes")

# results <- fread("/scratch/ejy4bu/drosophila/gowinda/results/gowinda_A_gene.txt", header=FALSE)

results <- read.delim("/scratch/ejy4bu/drosophila/gowinda/results/gowinda_AB_gene.txt", header=FALSE, col.names=cols)

significant = results %>% filter(FDR < 0.05)


# Significant hits
results[FDR < 0.05]
results[p.value < 0.05]