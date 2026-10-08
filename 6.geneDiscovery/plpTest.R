### FROM ALAN

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


### canddiate genes
    ggplot(data=nlp, aes(x=maf_mel, y=maf_sim)) + geom_tile() + facet_grid(~call)


    nlp.tsp.gene <- nlp[!is.na(nLocales_poly_proj_round) & genmap_score==1 & busco=="Complete", 
                        list(.N), 
                        list(gene, col, call, cosmo=as.factor(nLocales_poly_proj_round>75))]
    nlp.tsp.gene.len <- nlp[,list(maxAApos=max(aa_pos)), list(gene)]
    nlp.tsp.gene <- merge(nlp.tsp.gene, nlp.tsp.gene.len, by="gene")

    ### version 1
        nlp.tsp.gene.w <- dcast(nlp.tsp.gene, gene+cosmo~col+call, value.var="N")
        nlp.tsp.gene.w[,mk:=(missense_variant_tsp/synonymous_variant_tsp)/(missense_variant_dmel/synonymous_variant_dmel)]
        ggplot(data=nlp.tsp.gene.w, aes(x=cosmo, y=log2(mk))) + geom_boxplot()

    ### version 2
        nlp.tsp.gene.w2 <- dcast(nlp.tsp.gene[cosmo==T], gene+maxAApos~col+call, value.var="N", fun.aggregate=sum)
        nlp.tsp.gene.w2[,mk:=(missense_variant_tsp/synonymous_variant_tsp)/(missense_variant_dmel/synonymous_variant_dmel)]
        nlp.tsp.gene.w2[,mk2:=(missense_variant_con/synonymous_variant_con)/(missense_variant_dmel/synonymous_variant_dmel)]

        ggplot(data=nlp.tsp.gene.w2, aes(log2(mk))) + geom_histogram()


    fisher.test(table(nlp.tsp.gene.w2$mk>1, nlp.tsp.gene.w2$mk2>1))

    nlp.tsp.gene.w2[log2(mk)>2 & mk!=Inf]
    nlp.tsp.gene.w2[mk>2][missense_variant_tsp>5]
    nlp.tsp.gene.w2[missense_variant_con>1]

    nlp[gene=="Sardh"][call=="tsp"][col%like%"missense"]

    ggplot(data=nlp[call=="tsp"], aes(x=maf_sim, y=maf_mel)) + geom_point()
    ggplot(data=nlp[call=="con"], aes(x=maf_sim, y=maf_mel)) + geom_point()

    fisher.test(table(
    nlp[call=="tsp"]$maf_sim >.05 & nlp[call=="tsp"]$maf_mel>.05,
    nlp[call=="tsp"]$col))

    nlp.tsp.gene[cosmo==T & call=="tsp"][col%like%"missense"][N>10]



### just the cosmo SNPs
    source("/Users/alanbergland/fbtr_to_fbgn.R")
    map <- read_fbtr_map("~/fbgn_fbtr_fbpp_expanded_fb_2026_03.tsv")   # read once
    nlp[, fbgn := fbtr_to_fbgn(transcript, map)$fbgn]
    nlp[, gene_symbol := fbtr_to_fbgn(transcript, map)$gene_symbol]
    nlp

    enz <- fread("/Users/alanbergland/Dmel_enzyme_data_fb_2026_03.tsv", skip=4)
    enz.ag <- enz[,list(enz=.N),gene_id]
    nlp <- merge(nlp, enz.ag, by.x="fbgn", by.y="gene_id", all.x=T)



    nlp.tsp.gene2 <- nlp[maf_sim >.25 & maf_mel>.25][enz>0][!is.na(nLocales_poly_proj_round) & genmap_score==1 & nLocales_poly_proj_round>=75, 
                        list(.N), 
                        list(fbgn,gene, col, call, cosmo=as.factor(nLocales_poly_proj_round>75))]
    nlp.tsp.gene.len <- nlp[,list(maxAApos=max(aa_pos)), list(gene)]
    nlp.tsp.gene2 <- merge(nlp.tsp.gene2, nlp.tsp.gene.len, by="gene")
    nlp.tsp.gene.w3 <- dcast(nlp.tsp.gene2[cosmo==T], fbgn+gene+maxAApos~col+call, value.var="N", fun.aggregate=sum)
    candidates <- nlp.tsp.gene.w3[missense_variant_tsp>=1][,c("fbgn", "gene", "maxAApos", "missense_variant_tsp", "synonymous_variant_tsp")]$fbgn
    universe <- nlp.tsp.gene.w3[missense_variant_tsp>=0][,  c("fbgn", "gene", "maxAApos", "missense_variant_tsp", "synonymous_variant_tsp")]$fbgn
    length(candidates)
    length(universe)

    write.csv(candidates, "~/candidats.csv", quote=F, row.names=F)
    write.csv(universe, "~/universe.csv", quote=F, row.names=F)



### FILTERING: From ALAN
# gives: CG6415, CG4573, CG31759
    # as nuclear encoded mito genes
nlp.tsp.gene2 <- nlp[maf_sim >.25 & maf_mel>.25][enz>0][!is.na(nLocales_poly_proj_round) & genmap_score==1 & nLocales_poly_proj_round>=75 & busco=="Complete",
                        list(.N),
                        list(fbgn,gene, col, call, cosmo=as.factor(nLocales_poly_proj_round>75))]
    nlp.tsp.gene.len <- nlp[,list(maxAApos=max(aa_pos)), list(gene)]
    nlp.tsp.gene2 <- merge(nlp.tsp.gene2, nlp.tsp.gene.len, by="gene")
    nlp.tsp.gene.w3 <- dcast(nlp.tsp.gene2[cosmo==T], fbgn+gene+maxAApos~col+call, value.var="N", fun.aggregate=sum)
    candidates <- nlp.tsp.gene.w3[missense_variant_tsp>=1][,c("fbgn", "gene", "maxAApos", "missense_variant_tsp", "synonymous_variant_tsp")]$fbgn
    universe <- nlp.tsp.gene.w3[missense_variant_tsp>=0][,  c("fbgn", "gene", "maxAApos", "missense_variant_tsp", "synonymous_variant_tsp")]$fbgn
    length(candidates)
    length(universe)