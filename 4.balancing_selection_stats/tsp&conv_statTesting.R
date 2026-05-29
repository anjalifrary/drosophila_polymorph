
# anova_dt <- shared_dt[!is.na(PostMedian)]
# anova_dt[, classification := factor(classification)]

# model <- lm(PostMedian ~ classification, data = anova_dt)
# anova(model)



library(data.table)
library(ggplot2)
library(foreach)


nlp <- load("/project/berglandlab/multispecies_endemism/data/collectiveAnalysis_version3/Drosophila_melanogaster.10_03_2026.nlpTable.paramask.genmap.busco.Rdata")

tsp <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed_geva.rds")

# merge tsp file on chr and pos, include classification column
nlp <- merge(nlp, tsp[,c("chr", "pos", "classification")], by=c("chr", "pos"), all.x=T)

# geva ages 
age <- fread("/scratch/ejy4bu/drosophila/AlleleAges.VA.cm_GEVA.txt")