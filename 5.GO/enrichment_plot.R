library(data.table)
library(ggplot2)
library(foreach)
library(dplyr)


cols <- c("GO.id", "SimulatedGenes", "ObservedGenes", 
  "p.value", "FDR", 
  "N_GenesFound", "N_GenesMax", "N_GenesTotal",
  "Description", "Genes")

# results <- fread("/scratch/ejy4bu/drosophila/gowinda/results/gowinda_A_gene.txt", header=FALSE)

# results <- read.delim("/scratch/ejy4bu/drosophila/gowinda/results/gowinda_AB_gene.txt", header=FALSE, col.names=cols)
dir <- "/scratch/ejy4bu/drosophila/gowinda/MAF5/"

library(ontologyIndex)
go_file <- get_ontology("/scratch/ejy4bu/drosophila/gowinda/go.obo")

results_AB <- read.delim(paste0(dir, "gowinda_AB_snp_allBackground.txt"), header = FALSE, col.names = cols)
results_FGOPXY <- read.delim(paste0(dir, "gowinda_FGOPXY_snp_allBackground.txt"), header = FALSE, col.names = cols)

# add english description
results_AB$Description <- go_file$name[results_AB$GO.id]
results_FGOPXY$Description <- go_file$name[results_FGOPXY$GO.id]

# filter for FDR < 0.05
results_AB <- results_AB[results_AB$FDR < 0.05,]
results_FGOPXY <- results_FGOPXY[results_FGOPXY$FDR < 0.05,]

# add ontology - most interested in BP 
library(AnnotationDbi)
library(GO.db)
go_dt <- data.table(GO.id = keys(GO.db))
go_dt[, ontology := Ontology(GO.id)]
# gowinda_results from line 22, modify for each results file

results_AB <- merge(results_AB, go_dt, by="GO.id", all.x=T)
results_FGOPXY <- merge(results_FGOPXY, go_dt, by="GO.id", all.x=T)

results_AB <- results_AB[results_AB$ontology=="BP", ]
results_FGOPXY <- results_FGOPXY[results_FGOPXY$ontology=="BP", ]


# calculate enrichment and -log10FDR
# results_AB[, enrichment := (ObservedGenes - ExpectedGenes) / ExpectedGenes]
results_AB$enrichment <- results_AB$ObservedGenes / results_AB$SimulatedGenes
# results_AB[, log2_enrichment := log2(ObservedGenes / SimulatedGenes)]
results_AB$neglogFDR <- -log10(results_AB$FDR)

# results_FGOPXY[, enrichment := (ObservedGenes - ExpectedGenes) / ExpectedGenes]
results_FGOPXY$enrichment <- results_FGOPXY$ObservedGenes / results_FGOPXY$SimulatedGenes
# results_FGOPXY[, log2_enrichment := log2(ObservedGenes / SimulatedGenes)]
results_FGOPXY$neglogFDR <- -log10(results_FGOPXY$FDR)


results <- results_AB
# results <- results_FGOPXY

ggplot(results, aes(x = enrichment, y = neglogFDR)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = 1, linetype = 2) +
  theme_bw() +
  xlab("Observed / Simulated Genes") +
  ylab("-log10(FDR)") 
#   + 
#   geom_hline(yintercept = -log10(0.05), linetype = 2)



ggplot(results, aes(x = enrichment, y = neglogFDR)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("grey70", "red")) +
  theme_bw()


# significant line cut off (horizontal = sig cutoff, vert = no enrichment)
geom_hline(yintercept = -log10(0.05), linetype = 2)

# right = overrepresented GO 
# left = underrep
# top = statistically robust