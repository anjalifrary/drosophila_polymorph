library(data.table)

results <- fread("/scratch/ejy4bu/drosophila/gowinda/results/gowinda_A_gene.txt", header=FALSE)
setnames(results, c("GO.id", "SimulatedGenes", "ObservedGenes", 
  "p.value", "FDR", 
  "N_GenesFound", "N_GenesMax", "N_GenesTotal",
  "Description", "Genes"))


# Significant hits
results[FDR < 0.05]
results[p.value < 0.05]