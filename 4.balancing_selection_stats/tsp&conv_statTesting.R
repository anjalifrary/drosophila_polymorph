
# anova_dt <- shared_dt[!is.na(PostMedian)]
# anova_dt[, classification := factor(classification)]

# model <- lm(PostMedian ~ classification, data = anova_dt)
# anova(model)

### old - see plots_mean_nlp-xtx-age.R

library(data.table)
library(ggplot2)
library(foreach)


load("/project/berglandlab/multispecies_endemism/data/collectiveAnalysis_version3/Drosophila_melanogaster.10_03_2026.nlpTable.paramask.genmap.busco.Rdata")

tsp <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed_geva.rds")

# merge tsp file on chr and pos, include classification column
nlp <- merge(nlp, tsp[,c("chr", "pos", "classification")], by=c("chr", "pos"), all.x=T)

# geva ages 
age <- fread("/scratch/ejy4bu/drosophila/AlleleAges.VA.cm_GEVA.txt")

# merge nlp with age on chr and pos 
age[,chr:=tstrsplit(id, "\\.")[[1]]]
age[,pos:=position]
nlp <- merge(nlp, age, by=c("chr", "pos"), all.x=T)

nlp[,status:="not_tsp"]                                 # default status is not_tsp
nlp[classification%in%c("A", "B"), status:="tsp"]

# # ??? 
# nlp[is.na(classification),status2:="not_tsp"]
# nlp[classification%in%c("A", "B"), status2:="tsp"]

tsp_conv <- nlp[(classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y"))]

anova(lm(PostMean ~ classification, data = tsp_conv[!is.na(PostMean)]))

# iterated through different pairs (AB, FG, OP, XY) - all statistically significant
subset <- tsp_conv[(classification%in%c("A","B", "X", "Y"))]
subset[, class2 := classification]
subset[classification%in%c("F", "O", "X"), class2 := "A_pair"]
subset[classification%in%c("G", "P", "Y"), class2 := "B_pair"]

anova(lm(PostMean ~ class2, data = subset[!is.na(PostMean)]))


#####################################################################################
### graph mean(PostMean) vs class , with stderr bars 

plot_dt <- tsp_conv[!is.na(PostMean) & classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y"), 
    .(mean_age = mean(PostMean), 
    stderr_age = sd(PostMean)/sqrt(.N),
    n = .N)
    , by = classification]

ggplot(data = plot_dt, aes(x=classification, y=mean_age)) + geom_point(size=3) +
    geom_errorbar(aes(ymin = mean_age - stderr_age, ymax = mean_age + stderr_age), width = 0.2)

