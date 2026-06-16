library(data.table)
library(ggplot2)
library(foreach)


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
