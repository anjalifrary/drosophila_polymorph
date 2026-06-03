### scp aob2x@login.hpc.virginia.edu:/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed.rds ~/.

### libraries
    library(data.table)
    library(ggplot2)
    library(foreach)

    

### load data
    # load("/Users/alanbergland/Downloads/Drosophila_melanogaster.10_03_2026.nlpTable.paramask.genmap.busco.Rdata") ###   save(nlp.mel, file="/project/berglandlab/multispecies_endemism/data/Dmel/Dmelanogaster_nlp.Rdata")
    # load("/scratch/ejy4bu/drosophila/Drosophila_melanogaster.10_03_2026.nlpTable.Rdata")
    load("/project/berglandlab/multispecies_endemism/data/collectiveAnalysis_version3/Drosophila_melanogaster.10_03_2026.nlpTable.paramask.genmap.busco.Rdata")
    tsp <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed_geva.rds")

    nlp <- merge(nlp, tsp[,c("chr", "pos", "classification")], by=c("chr", "pos"), all.x=T)

    age <- fread("/scratch/ejy4bu/drosophila/AlleleAges.VA.cm_GEVA.txt")
    age[,chr:=tstrsplit(id, "\\.")[[1]]]
    age[,pos:=position]
    nlp <- merge(nlp, age, by=c("chr", "pos"), all.x=T)

    nlp[,status:="not_tsp"]
    nlp[classification%in%c("A", "B"), status:="tsp"]

    nlp[is.na(classification),status2:="not_tsp"]
    nlp[classification%in%c("A", "B"), status2:="tsp"]


    # Alan's plots
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


# summaries
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


###########################################################################################
# fisher's exact test on A/B vs the natural pairs
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

# AB vs FG
count <- c(classA_count, classF_count, classB_count, classG_count)
ft <- fisher.test(matrix(count, nrow=2, byrow=T))  # odds ratio = 5.582855 ; p-value < 2.2e-16

# AB vs OP
count <- c(classA_count, classO_count, classB_count, classP_count)
ft <- fisher.test(matrix(count, nrow=2, byrow=T))
    # odds ratio =  0.1437368 ; p-value < 2.2e-16

# AB vs XY
count <- c(classA_count, classX_count, classB_count, classY_count)
ft <- fisher.test(matrix(count, nrow=2, byrow=T))
    # odds ratio = 9.156533 ; p-value < 2.2e-16

# AB vs FO, GP
count <- c(classA_count, sum(classF_count+classO_count), classB_count, sum(classG_count+classP_count))
ft <- fisher.test(matrix(count, nrow=2, byrow=T))
    # odds ratio = 3.307256 ; p-value < 2.2e-16

# AB vs FOX, GPY
count <- c(classA_count, sum(classF_count,classO_count,classX_count), classB_count, sum(classG_count, classP_count, classY_count))
ft <- fisher.test(matrix(count, nrow=2, byrow=T))
    # odds ratio = 3.521055 ; p-value < 2.2e-16


### bar graph of relative proportion of shared vs independent (colored by NS/Syn)
df <- data.frame(
    Category = c("Shared", "Independent", "Shared", "Independent"),
    Type = c("NS", "NS", "Syn", "Syn"),
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



# from Alan, 5/27
### libraries
    library(data.table)
    library(ggplot2)
    library(foreach)
    library(doMC)
    registerDoMC(10)

### load data
    load("/Users/alanbergland/Downloads/Drosophila_melanogaster.10_03_2026.nlpTable.paramask.genmap.busco.Rdata") ###   save(nlp.mel, file="/project/berglandlab/multispecies_endemism/data/Dmel/Dmelanogaster_nlp.Rdata")

    tsp <- readRDS("~/subset_qualVar_ofInterest_classed.rds")
 
    nlp <- merge(nlp, tsp[,c("chr", "pos", "classification")], by=c("chr", "pos"), all.x=T)

    age <- fread("/Users/alanbergland/Downloads/AlleleAges.VA.cm_GEVA.txt")
    age[,chr:=tstrsplit(id, "\\.")[[1]]]
    age[,pos:=position] 
    nlp <- merge(nlp, age, by=c("chr", "pos"), all.x=T)

    nlp[,status:="not_tsp"]
    nlp[classification%in%c("A", "B"), status:="tsp"]
    
    nlp[is.na(classification),status2:="not_tsp"]
    nlp[classification%in%c("A", "B"), status2:="tsp"]

### age clustering
    library(mclust)
    fit_mclust <- Mclust(nlp$PostMean, G=2)    
    summary(fit_mclust, parameters = TRUE)
    nlp[!is.na(PostMean), ageClust:=fit_mclust$classification]
    ggplot(data=nlp, aes(PostMean, group=ageClust, fill=ageClust)) + geom_histogram()

### summarize

    nlp.ag <- nlp[genmap_score==1,
                list(prTSP=mean(status=="tsp", na.rm=T), tspN=sum(!is.na(status)), 
                    prOld=mean(ageClust==2, na.rm=T), oldAge=mean(PostMean[ageClust==2], na.rm=T), 
                    age=mean(PostMean, na.rm=T), ageN=sum(!is.na(PostMean)), ageSD=sd(PostMean, na.rm=T)),
                list(nLocales_poly=nLocales_poly, col=col)]
    
    nlp.ag.w <- dcast(nlp.ag, nLocales_poly~col, value.var=c("prTSP", "age"))
    nlp.ag.w[,rr:=(prTSP_missense_variant)/(prTSP_synonymous_variant)]
    nlp.ag.w[,ageDelta:=age_missense_variant - age_synonymous_variant]

    #ggplot(data=nlp.ag, aes(x=nLocales_poly, y=prTSP, group=col, color=col)) + geom_line()


### save
    save(nlp.ag, nlp.ag.w, file="/Users/alanbergland/Documents/GitHub/misc/1000G/done/Dmel/Dmel_figures/Figure4_balancing/parts/tsp_and_age.Rdata")



    ggplot(data=nlp.ag.w, aes(x=nLocales_poly, y=log2(rr))) + geom_line() +
    ylab("log2(Relative Rate):\nPr(TSP|NS)/Pr(TSP|S)") + theme_bw() +
    geom_smooth(method="lm") + ggtitle("Trans-specificity with D. simulans")
    
    ggplot(data=nlp.ag.w, aes(x=nLocales_poly, y=(ageDelta))) + geom_line()

    ggplot(data=nlp.ag, aes(x=nLocales_poly, y=prTSP, color=col, group=col)) + geom_line()
    ggplot(data=nlp.ag, aes(x=nLocales_poly, y=(age), color=col, group=col)) + geom_line()
    ggplot(data=nlp.ag, aes(x=nLocales_poly, y=prOld, color=col, group=col)) + geom_line()
    ggplot(data=nlp.ag, aes(x=nLocales_poly, y=oldAge, color=col, group=col)) + geom_line()




### for Anjali
    tsp.ag <- nlp[,list(prOld=mean(ageClust==2, na.rm=T), 
                       oldAge=mean(PostMean[ageClust==2], na.rm=T),
                       age=mean(PostMean, na.rm=T),
                       nlp=mean(nLocales_poly, na.rm=T), 
                       nlpOld=mean(nLocales_poly[ageClust==2], na.rm=T)), 
                    list(classification, col)]

    tsp.ag.
    ggplot(data=tsp.ag[classification%in%c("A", "F", "O", "B", "G", "P")], aes(x=classification, y=prOld, color=col)) + geom_point() + facet_grid(~col, scales="free_x")
    ggplot(data=tsp.ag[classification%in%c("A", "F", "O", "B", "G", "P")], aes(x=classification, y=oldAge, color=col)) + geom_point() + facet_grid(~col, scales="free_x")
    ggplot(data=tsp.ag[classification%in%c("A", "F", "O", "B", "G", "P")], aes(x=classification, y=age, color=col)) + geom_point() + facet_grid(~col, scales="free_x")
    ggplot(data=tsp.ag[classification%in%c("A", "F", "O", "B", "G", "P")], aes(x=classification, y=nlp, color=col)) + geom_point() + facet_grid(~col, scales="free_x")
    ggplot(data=tsp.ag[classification%in%c("A", "F", "O", "B", "G", "P")], aes(x=classification, y=nlpOld, color=col)) + geom_point() + facet_grid(~col, scales="free_x")


    summary(lm(PostMean~classification, nlp[classification%in%c("A", "F", "O")]))
    summary(glm(I(ageClust==2)~classification, nlp[classification%in%c("A", "F", "O")], family=binomial()))


### gene level OR
    buffer <- 0
    nlp.ag <- nlp[genmap_score==1 ,list(tsp_NS=sum(classification=="A" & global_af>.01, na.rm=T)+buffer, 
              tsp_SYN=sum(classification=="B" & global_af>.01, na.rm=T)+buffer,
              indp_NS=sum(classification%in%c("F", "O"), na.rm=T)+buffer, 
              indp_SYN=sum(classification%in%c("G", "P"), na.rm=T)+buffer, 
              poly_NS=sum(col%like%"missense" & status=="not_tsp" & global_af>.01, na.rm=T)+buffer,
              poly_SYN=sum(col%like%"synonymous" & status=="not_tsp" & global_af>.01, na.rm=T)+buffer,
              age=mean(PostMean, na.rm=T),
              prOld=mean(ageClust==2, na.rm=T),
              nlp=mean(nLocales_poly, na.rm=T)),
        list(gene)]

    nlp.ag[,or_tsp:=(tsp_NS/tsp_SYN)/(poly_NS/poly_SYN)]
    nlp.ag[,or_indp:=(indp_NS/indp_SYN)/(poly_NS/poly_SYN)]

    nlp.ag[or_tsp>2][or_tsp!=Inf][(tsp_NS+tsp_SYN)>10][order(or_tsp)]
    nlp.ag[or_indp>3][or_indp!=Inf][(indp_NS+indp_SYN)>10][order(or_indp)]

    hist(log2(nlp.ag[or_tsp!=Inf & or_tsp!=0]$or_tsp))

    hist(log2(nlp.ag$or))

    ggplot(data=nlp.ag, aes(x=log2(or_tsp), y=log2(or_indp))) + geom_point()
    ggplot(data=nlp.ag, aes(x=nlp, y=log2(or_indp))) + geom_point(alpha=.5)
    ggplot(data=nlp.ag, aes(x=prOld, y=log2(or_tsp))) + geom_point()

    summary(lm(log2(or_indp)~nlp, nlp.ag))
    summary(lm(log2(or_tsp)~log2(or_indp), nlp.ag[or_tsp!=0 & or_tsp!=Inf & or_indp!=0 & or_indp!=Inf]))nlp

    nlp.ag[or_indp==Inf][indp_NS>5]


    tab <- table(nlp.ag$or_tsp>2, nlp.ag$or_indp>2)
    fisher.test(tab)
    nlp.ag[gene%like%"Adh"]

    ggplot(data=tsp.ag[classification%in%c("B", "G", "P")], aes(x=classification, y=prOld, color=col)) + geom_point()

    ggplot(data=tsp.ag, aes(x=classification, y=oldAge, color=col)) + geom_point()












scp aob2x@login.hpc.virginia.edu:/project/berglandlab/anjali/drosophila_polymorphism/classification/classification_table_final.csv ~/.

library(data.table)
library(ggplot2)
options(width = 170)
dat <- fread("~/classification_table_final.csv")
dat











### load windows of interest
  load("/Users/alanbergland/Documents/GitHub/misc/1000G/done/Dmel/Dmel_figures/Figure4_balancing/gwas/Genomic_windows_EdemismOutliers_Clusters_and_BGD.save.Rdata")
  sw <- Genomic_windows_EdemismOutliers_Clusters_and_BGD.save
  sw <- as.data.table(sw)
  swl <- foreach(j=1:dim(sw)[1], .combine="rbind")%dopar%{
    # j <- 1
    data.table(winId=paste(sw[j]$SigWinId, sw[j]$cluster, sep="_"), chr=sw[j]$chr, pos=sw[j]$begin:sw[j]$end)
  }
  setkey(swl, chr, pos)

    age <- fread("/Users/alanbergland/Downloads/AlleleAges.VA.cm_GEVA.txt")
    age[,chr:=tstrsplit(id, "\\.")[[1]]]
    age[,pos:=position] 

   library(mclust)
    fit_mclust <- Mclust(age$PostMean, G=2)    
    summary(fit_mclust, parameters = TRUE)
    age[!is.na(PostMean), ageClust:=fit_mclust$classification]
   age_hist <- ggplot(data=age, aes(PostMean, group=ageClust, fill=as.factor(ageClust))) + geom_histogram( alpha=.5, bins=100)



setkey(age, chr, pos)
agew <- merge(age, swl)
agew.ag <- agew[,list(age=mean(PostMean, na.rm=T), age.sd=sd(PostMean, na.rm=T), prOld=mean(ageClust==2), oldAge=mean(PostMean[ageClust==2], na.rm=T), oldAge.sd=sd(PostMean[ageClust==2], na.rm=T)), list(winId)]


agew.ag[,class:=tstrsplit(winId, "_")[[4]]]
o.ag <- agew.ag[,list(age=mean(age), age.sd=sd(age), .N, prOld=mean(prOld), prOld.sd=sd(prOld), oldAge=mean(oldAge), oldAge.sd=sd(oldAge)), list(class)]
o.ag[,se:=age.sd/sqrt(N)]
o.ag[,prOldSE:=prOld.sd/sqrt(N)]
o.ag[,oldAgeSE:=oldAge.sd/sqrt(N)]


age_plot <- ggplot(data=agew.ag, aes(x=class, y=age)) + geom_beeswarm() + geom_point(data=o.ag, aes(x=class, y=age), color="red", size=4) + geom_linerange(data=o.ag, aes(ymin=age-2*se, ymax=age+2*se), color="red") + theme_bw()
prOld_plot <- ggplot(data=agew.ag, aes(x=class, y=prOld)) + geom_beeswarm() + geom_point(data=o.ag, aes(x=class, y=prOld), color="red", size=4) + geom_linerange(data=o.ag, aes(ymin=prOld-2*prOldSE, ymax=prOld+2*prOldSE), color="red") + theme_bw()
oldAge_plot <- ggplot(data=agew.ag, aes(x=class, y=oldAge)) + geom_beeswarm() + geom_point(data=o.ag, aes(x=class, y=oldAge), color="red", size=4) + geom_linerange(data=o.ag, aes(ymin=oldAge-2*oldAgeSE, ymax=oldAge+2*oldAgeSE), color="red") + theme_bw()


age_plot <- ggplot(data=o.ag, aes(x=class, y=age)) + geom_linerange(aes(ymin=age-2*se, ymax=age+2*se)) + geom_point() + theme_bw()
prOld_plot <- ggplot(data=o.ag, aes(x=class, y=prOld)) + geom_linerange(aes(ymin=prOld-2*prOldSE, ymax=prOld+2*prOldSE)) + geom_point() + theme_bw()
oldAge_plot <- ggplot(data=o.ag, aes(x=class, y=oldAge)) + geom_linerange(aes(ymin=oldAge-2*oldAgeSE, ymax=oldAge+2*oldAgeSE)) + geom_point() + theme_bw()

age_hist + age_plot + prOld_plot + oldAge_plot
