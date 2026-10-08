# 10/8/26:

### libraries
    library(data.table)
    library(ggplot2)
    library(patchwork)

### load data
    load("/Users/alanbergland/Drosophila_melanogaster.06_10_2026.nlpTable.paramask.genmap.busco.repeatmask.wmdust.projection.Rdata")

    nlp.ag <- nlp[,
                  list(nNS=sum(col%like%"missense"), nS=sum(col%like%"synonymous"),
                        lNS=sum(nonsynonymous), lS=sum(synonymous)),
                  list(nLocales_poly_proj_round=nLocales_poly_proj_binomial)]

    nlp.ag[,propNS:=nNS/(nNS+nS)]
    nlp.ag[,N:=(nNS+nS)]
    nlp.ag[,pnps:=(nNS/lNS)/(nS/lS)]

    ggplot(data=nlp.ag, aes(x=nLocales_poly_proj_round, y=propNS)) + geom_line()
    ggplot(data=nlp.ag, aes(x=nLocales_poly_proj_round, y=pnps)) + geom_line()

    ggplot(data=nlp.ag.w, aes(x=nLocales_poly_proj_round, y=log10(N))) + geom_line()


### tsps
    tsp <- readRDS("/Users/alanbergland/dsim3.signor.DGRP2.source_BCM-HGSC.shared.bothMelSim.classed.rds")
    str(tsp)

    tsp.small <- tsp[, c("classification", "maf_sim", "maf_mel", "chr_dm6", "pos_dm6")]
    setnames(tsp.small, c("chr_dm6", "pos_dm6"), c("chr", "pos"))
    setkey(tsp.small, chr, pos)
    setkey(nlp, chr, pos)
    nlp <- merge(nlp, tsp.small, all.x=T)
    nlp[classification%in%c("A","B"),call:="tsp"]
    nlp[classification%in%c("F", "G", "O", "P", "X", "Y"), call:="con"]
    nlp[is.na(call), call:="dmel"]

    nlp.tsp <- nlp[, list(nNS=sum(col%like%"missense"), nS=sum(col%like%"synonymous"),
                          lNS=sum(nonsynonymous), lS=sum(synonymous)),
                     list(nLocales_poly_proj_round=nLocales_poly_proj_round,
                          call=call)]
    nlp.tsp[,N:=nNS+nS]
    nlp.tsp[,propNS:=nNS/(nNS+nS)]
    nlp.tsp[,pnps:=(nNS/lNS)/(nS/lS)]

    nlp.tsp[,fracSNPs:=N/sum(nlp.tsp$N)]

    nlp.tsp.w <- dcast(nlp.tsp, nLocales_poly_proj_round~call, value.var="N")
    nlp.tsp.w[,tsp_v_dmel:=log2(tsp/dmel)]
    nlp.tsp.w[,con_v_dmel:=log2(con/dmel)]
    nlp.tsp.rel <- melt(nlp.tsp.w, id.vars="nLocales_poly_proj_round", measure.vars=c("tsp_v_dmel", "con_v_dmel"))


    nlp.tsp.2w <- dcast(nlp.tsp, nLocales_poly_proj_round~call, value.var="pnps")
    nlp.tsp.2w[,tsp_v_dmel:=log2(tsp/dmel)]
    nlp.tsp.2w[,con_v_dmel:=log2(con/dmel)]
    nlp.tsp.2rel <- melt(nlp.tsp.2w, id.vars="nLocales_poly_proj_round", measure.vars=c("tsp_v_dmel", "con_v_dmel"))



    a <- ggplot(data=nlp.tsp, aes(x=nLocales_poly_proj_round, y=log10(N), group=call, color=call)) + geom_line()
    b <- ggplot(data=nlp.tsp.rel, aes(x=nLocales_poly_proj_round, y=value, group=variable, color=variable)) + geom_line() + ylab("log2(# class/# dmel)")
    c <- ggplot(data=nlp.tsp, aes(x=nLocales_poly_proj_round, y=pnps, group=call, color=call)) + geom_line()
    d <- ggplot(data=nlp.tsp.2rel, aes(x=nLocales_poly_proj_round, y=value, group=variable, color=variable)) + geom_line() + ylab("log2(pnps class/pnps dmel)")

    layout <-"
    AB
    CD"
    a + b + c + d + plot_layout(design=layout)