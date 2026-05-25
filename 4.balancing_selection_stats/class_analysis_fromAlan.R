### scp aob2x@login.hpc.virginia.edu:/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed.rds ~/.

### libraries
    library(data.table)
    library(ggplot2)
    library(foreach)

### load data
    # load("/Users/alanbergland/Downloads/Drosophila_melanogaster.10_03_2026.nlpTable.paramask.genmap.busco.Rdata") ###   save(nlp.mel, file="/project/berglandlab/multispecies_endemism/data/Dmel/Dmelanogaster_nlp.Rdata")
    load("/scratch/ejy4bu/drosophila/Drosophila_melanogaster.10_03_2026.nlpTable.Rdata")
    tsp <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed_geva.rds")

    nlp <- merge(nlp, tsp[,c("chr", "pos", "classification")], by=c("chr", "pos"), all.x=T)

    age <- fread("/scratch/ejy4bu/drosophila/AlleleAges.VA.cm_GEVA.txt")
    age[,chr:=tstrsplit(id, "\\.")[[1]]]
    age[,pos:=position]
    nlp <- merge(nlp, age, by=c("chr", "pos"), all.x=T)

    nlp[,status:="not_tsp"]
    nlp[classification%in%c("A", "B"), status:="tsp"]
    ggplot(data=nlp, aes(x=status, y=nLocales_poly)) + geom_boxplot()
    # ggplot(data=nlp, aes(x=status, y=PostMean)) + geom_boxplot()
    ggplot(data=nlp, aes(x=status, y=PostMode)) + geom_boxplot()
    # ggplot(data=nlp, aes(x=nLocales_poly, y=PostMean, color=status)) + geom_point()
    ggplot(data=nlp, aes(x=nLocales_poly, y=PostMode, color=status)) + geom_point()
    # ggplot(data=nlp, aes(x=classification, y=nLocales_poly)) + geom_boxplot()
    ggplot(data=nlp[!is.na(classification)], aes(x=classification, y=nLocales_poly)) + geom_boxplot()


    nlp[,list(pr=mean(nLocales_poly>100), .N), list(col, classification)][order(-pr)]


    fisher.test(table(nlp$status, nlp$nLocales_poly>100))           # odds ratio = 3.78

    fisher.test(table(nlp[nLocales_poly>100 & PostMean<5e5]$status, # odds ratio = 2.237165
                    nlp[nLocales_poly>100 & PostMean<5e5]$col))
    fisher.test(table(nlp[nLocales_poly>100 & PostMode<5e5]$status,
                    nlp[nLocales_poly>100 & PostMode<5e5]$col))     # odds ratio = 2.25183

    fisher.test(table(nlp[nLocales_poly>100 & PostMean>2e6]$status, # odds ratio = 2.762643 
                    nlp[nLocales_poly>100 & PostMean>2e6]$col))
    fisher.test(table(nlp[nLocales_poly>100 & PostMode>2e6]$status,
                    nlp[nLocales_poly>100 & PostMode>2e6]$col))     # odds ratio = 2.764495 

    nlp[nLocales_poly>100 & classification=="F" & col%like%"missense"]
    nlp[nLocales_poly>100 & classification=="O" & col%like%"missense"]

    summary(glm(I(nLocales_poly>100)~classification, data=nlp, family=binomial()))


# ###
    nlp.ag <- nlp[,
                 list(prTSP=mean(status=="tsp", na.rm=T), N=sum(status=="tsp", na.rm=T), age=mean(PostMean, na.rm=T)),
                 list(nLocales_poly=nLocales_poly, col=col)]

    # nlp.ag <- nlp[,
    #              list(prTSP=mean(status=="tsp", na.rm=T), N=sum(status=="tsp", na.rm=T), age=mean(PostMode, na.rm=T)),
    #              list(nLocales_poly=nLocales_poly, col=col)]

    nlp.ag.w <- dcast(nlp.ag, nLocales_poly~col, value.var="prTSP")
    nlp.ag.w[,rr:=(missense_variant)/(synonymous_variant)] # relative rate 
    #ggplot(data=nlp.ag, aes(x=nLocales_poly, y=prTSP, group=col, color=col)) + geom_line()

    tsp_fig <- ggplot(data=nlp.ag.w, aes(x=nLocales_poly, y=log2(rr))) + geom_line() +
    ylab("log2(Relative Rate):\nPr(TSP|NS)/Pr(TSP|S)") + theme_bw() +
    geom_smooth(method="lm") + ggtitle("Trans-specificity with D. simulans")


# Try reconstructing this type of analysis for these sets of pairs:
# `fisher.test(matrix(c(13906, 66171, 9323, 28508), nrow=2, byrow=T))

###########################################################################################
###########################################################################################

class_table <- read.csv("/project/berglandlab/anjali/drosophila_polymorphism/classification/classification_table_final.csv")

classA_count <- class_table[class_table$Classification=="A", "Count"]
classB_count <- class_table[class_table$Classification=="B", "Count"]

classF_count <- class_table[class_table$Classification=="F", "Count"]
classG_count <- class_table[class_table$Classification=="G", "Count"]

classO_count <- class_table[class_table$Classification=="O", "Count"]
classP_count <- class_table[class_table$Classification=="P", "Count"]

classX_count <- class_table[class_table$Classification=="X", "Count"]
classY_count <- class_table[class_table$Classification=="Y", "Count"]

nonsyn <- list(classA_count, classF_count, classO_count, classX_count)
syn <- list(classB_count, classG_count, classP_count, classY_count)

test1 <- matrix(c(classA_count, classF_count, classB_count, classG_count), nrow=2, byrow=T)
    #       (shared)   (ind)
    # (NS)  13906 (A)  1073  (F)
    # (syn) 66171 (B)  28508 (G)

count <- c(classA_count, classF_count, classB_count, classG_count)
ft <- fisher.test(matrix(count, nrow=2, byrow=T))  # odds ratio = 5.582855 ; p-value < 2.2e-16

count <- c(classA_count, classO_count, classB_count, classP_count)
ft <- fisher.test(matrix(count, nrow=2, byrow=T))
    # odds ratio =  0.1437368 ; p-value < 2.2e-16

count <- c(classA_count, classX_count, classB_count, classY_count)
ft <- fisher.test(matrix(count, nrow=2, byrow=T))
    # odds ratio = 9.156533 ; p-value < 2.2e-16

count <- c(classA_count, sum(classF_count+classO_count), classB_count, sum(classG_count+classP_count))
ft <- fisher.test(matrix(count, nrow=2, byrow=T))
    # odds ratio = 3.307256 ; p-value < 2.2e-16

count <- c(classA_count, sum(classF_count,classO_count,classX_count), classB_count, sum(classG_count, classP_count, classY_count))
ft <- fisher.test(matrix(count, nrow=2, byrow=T))
    # odds ratio = 3.521055 ; p-value < 2.2e-16


### bar graph of relative proportion of shared vs independent (colored by NS/Syn)
df <- data.frame(
    Category = c("Shared", "Shared", "Independent", "Independent"),
    Type = c("NS", "Syn", "NS", "Syn"),
    Count = count
)
ggplot(df, aes(x=Category, y=Count, fill=Type)) +
  geom_bar(stat="identity", position="fill") +
  ylab("Proportion") +
  theme_bw()

### odds ratio graph with error bar 
plot_df <- data.frame(
  comparison = "Shared vs Independent",
  odds_ratio = ft$estimate,
  lower = ft$conf.int[1],
  upper = ft$conf.int[2]
)
ggplot(data=plot_df, aes(x=comparison, y = odds_ratio)) + geom_point(size=4) + geom_hline(yintercept=1) + 
    ylab("Odds Ratio") + geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)
