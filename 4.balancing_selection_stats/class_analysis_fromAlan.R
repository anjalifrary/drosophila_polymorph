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
    ggplot(data=nlp, aes(x=status, y=PostMean)) + geom_boxplot()
    ggplot(data=nlp, aes(x=nLocales_poly, y=PostMean, color=status)) + geom_point()
    ggplot(data=nlp, aes(x=classification, y=nLocales_poly)) + geom_boxplot()
    ggplot(data=nlp[!is.na(classification)], aes(x=classification, y=nLocales_poly)) + geom_boxplot()

    nlp[,list(pr=mean(nLocales_poly>100), .N), list(col, classification)][order(pr)]


    fisher.test(table(nlp$status, nlp$nLocales_poly>100))

    fisher.test(table(nlp[nLocales_poly>100 & PostMean<5e5]$status,
                    nlp[nLocales_poly>100 & PostMean<5e5]$col))

    fisher.test(table(nlp[nLocales_poly>100 & PostMean>2e6]$status,
                    nlp[nLocales_poly>100 & PostMean>2e6]$col))

    nlp[nLocales_poly>100 & classification=="K" & col%like%"missense"]
    nlp[nLocales_poly>100 & classification=="K" & col%like%"missense"]

    summary(glm(I(nLocales_poly>100)~classification, data=nlp, family=binomial()))


###
    nlp.ag <- nlp[genmap_score==1,
                 list(prTSP=mean(status=="tsp", na.rm=T), N=sum(status=="tsp", na.rm=T), age=mean(PostMean, na.rm=T)),
                 list(nLocales_poly=nLocales_poly, col=col)]

    nlp.ag.w <- dcast(nlp.ag, nLocales_poly~col, value.var="prTSP")
    nlp.ag.w[,rr:=(missense_variant)/(synonymous_variant)]
    #ggplot(data=nlp.ag, aes(x=nLocales_poly, y=prTSP, group=col, color=col)) + geom_line()

    tsp_fig <- ggplot(data=nlp.ag.w, aes(x=nLocales_poly, y=log2(rr))) + geom_line() +
    ylab("log2(Relative Rate):\nPr(TSP|NS)/Pr(TSP|S)") + theme_bw() +
    geom_smooth(method="lm") + ggtitle("Trans-specificity with D. simulans")


# Try reconstructing this type of analysis for these sets of pairs:
# `fisher.test(matrix(c(13906, 66171, 9323, 28508), nrow=2, byrow=T))

