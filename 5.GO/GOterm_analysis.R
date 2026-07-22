library(data.table)
library(ggplot2)
library(foreach)


### GO Sem sim 
# https://yulab-smu.top/biomedical-knowledge-mining-book/011-GOSemSim.html
## citation
# Guangchuang Yu, Fei Li, Yide Qin, Xiaochen Bo, Yibo Wu and Shengqi Wang.
# GOSemSim: an R package for measuring semantic similarity among GO terms
# and gene products. Bioinformatics. 2010, 26(7):976-978

library(AnnotationDbi)
library(org.Dm.eg.db)
library(GOSemSim)
dmGO <- godata(annoDb = 'org.Dm.eg.db', ont="BP")
# User can set computeIC=FALSE if they only want to use Wang’s method.

# test similarity score between 0 and 1 [immune response vs immune system process]
# 1 = same
# 0 = completely different
goSim("GO:0006955", "GO:0002376", semData=dmGO, measure="Wang")
# gives 0.63

# get results tables
    library(dplyr)


    cols <- c("GO.id", "SimulatedGenes", "ObservedGenes", 
    "p.value", "FDR", 
    "N_GenesFound", "N_GenesMax", "N_GenesTotal",
    "Description", "Genes")

    # results <- fread("/scratch/ejy4bu/drosophila/GO/gowinda/results/gowinda_A_gene.txt", header=FALSE)

    # results <- read.delim("/scratch/ejy4bu/drosophila/GO/gowinda/results/gowinda_AB_gene.txt", header=FALSE, col.names=cols)
    dir <- "/scratch/ejy4bu/drosophila/GO/gowinda/MAF5/new_6-29-26/results/"

    library(ontologyIndex)
    go_file <- get_ontology("/scratch/ejy4bu/drosophila/GO/gowinda/go.obo")

    results_AB <- read.delim(paste0(dir, "gowinda_AB_gene_allBackground.txt"), header = FALSE, col.names = cols)
    results_FGOPXY <- read.delim(paste0(dir, "gowinda_FGOPXY_gene_allBackground.txt"), header = FALSE, col.names = cols)

    # add english description
    results_AB$Description <- go_file$name[results_AB$GO.id]
    results_FGOPXY$Description <- go_file$name[results_FGOPXY$GO.id]

    # add ontology - most interested in BP 
    library(GO.db)
    go_dt <- data.table(GO.id = keys(GO.db))
    go_dt[, ontology := Ontology(GO.id)]
    # gowinda_results from line 22, modify for each results file

    results_AB <- merge(results_AB, go_dt, by="GO.id", all.x=T)
    results_FGOPXY <- merge(results_FGOPXY, go_dt, by="GO.id", all.x=T)



AB_terms <- unique(results_AB$GO.id[results_AB$FDR<0.05])
FGOPXY_terms <- unique(results_FGOPXY$GO.id[results_FGOPXY$FDR<0.05])
length(AB_terms) # 1149
length(FGOPXY_terms) # 145

# mgoSim() calculates semantic similarity between two sets of GO terms 
AB_sim <- mgoSim(AB_terms, AB_terms, semData=dmGO, measure="Wang", combine=NULL)
FGOPXY_sim <- mgoSim(FGOPXY_terms, FGOPXY_terms, semData=dmGO, measure="Wang", combine=NULL)

all_sim <- mgoSim(AB_terms, FGOPXY_terms, semData=dmGO, measure="Wang", combine="BMA")
all_sim
# all_sim = 0.331

hc_AB <- hclust(as.dist(1 - AB_sim))
plot(hc_AB, cex = 0.4)

hc_FGOPXY <- hclust(as.dist(1-FGOPXY_sim))
plot(hc_FGOPXY, labels = FALSE)

hc_all <- hclust(as.dist())

### trying to view clusters
cl_AB <- cutree(hc_AB, k = 10)
cluster_AB <- data.table(
    GO.id = names(cl_AB),
    cluster = cl_AB
)

cluster_AB[
    results_AB[, .(GO.id, Description)],
    on = "GO.id"
][order(cluster)]

# mclusterSim() calculates semantic similarity between two sets of gene clusters
# clusterSim() between two gene clusters

# GO classification based on GO distribution at a specific level
library(clusterProfiler) # G Yu. Thirteen years of clusterProfiler. The Innovation. 2024,5(6):100722
library(org.Dm.eg.db)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
head(gene)


ggo <- groupGO(gene     = gene,
               OrgDb    = org.Dm.eg.db,
               keytype = "ENTREZID"
               ont      = "BP",
               level    = 3,
               readable = TRUE)

head(ggo)

### basically a repeat of gowinda but gene-level, not snp-level
# # GO over-represenation analysis
# ego <- enrichGO(gene          = gene,
#                 universe      = names(geneList),
#                 OrgDb         = org.Dm.eg.db,
#                 ont           = "BP",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.05,
#         readable      = TRUE)
# head(ego)

# generate var file from the geographic_clines.R file
var <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_ABFGOPXY_mafmel5_xtx_061626.rds")

carb_genes <- unique(unlist(strsplit(results_AB[GO.id=="GO:0005975", "Genes"][[1]], ",")))

var[, gene_id_mel := tolower(gene_id_mel)]

carb_var <- var[gene_id_mel%in%(carb_genes)]
# conv  tsp 
#  209  720 

#   A   B   F   G   O   P   X   Y 
#   79 641   7 173   6   6   1  16 

#   missense_variant synonymous_variant 
#   93                833 

# 266 cosmopolitan 
# 99 unique genes 

# filtered for missense and cosmopolitan (nLocales_poly > 100) = all class A
    # 20 variants, 17 genes



### get genes in GO ids like "metabolic"
metabolic_genes <- rbindlist(list(
    results_AB[ontology=="BP" & grepl("metabol", Description, ignore.case = TRUE), .(Genes)],
    results_FGOPXY[ontology=="BP" & grepl("metabol", Description, ignore.case = TRUE), .(Genes)])
)
metabolic_genes <- unlist(strsplit(
    metabolic_genes$Genes, ","
))

subset <- var[gene_id_mel%in%(metabolic_genes)]
subset <- subset[nLocales_poly>100 & class=="tsp"]
View(subset[col%like%c("missense")])

### immune genes
immune_genes <- rbindlist(list(
    results_AB[ontology=="BP" & grepl("immune", Description, ignore.case = TRUE), .(Genes)],
    results_FGOPXY[ontology=="BP" & grepl("immune", Description, ignore.case = TRUE), .(Genes)])
)

immune_genes <- unlist(strsplit(
    immune_genes$Genes, ","
))


subset <- var[gene_id_mel%in%(immune_genes)]
subset <- subset[poly_af>0.2]
