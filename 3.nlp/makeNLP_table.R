library(SeqArray)
library(data.table)

# meta <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/dest_v2.samps_24Aug2024.xa.csv")
meta <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/simulans_pooled.meta.use.csv")

samps <- meta[, .(
    sample = identifier, # sim
    pop = sub("_[^_]+_[0-9]{4}-[0-9]{2}-[0-9]{2}$", "", identifier),
    # sample = sampleId, # mel
    # pop = sub("_[^_]+_[0-9]{4}-[0-9]{2}-[0-9]{2}$", "", sampleId),
    nInd = nFlies,
    country,
    city,
    lat,
    long
)]

### common functions

fis <- function(nHet, nTot, p=.5) {
numerator <- (nHet/nTot)
denominator <- (nTot/(nTot-1))*(2*p*(1-p) - (nHet/nTot)/(2*nTot))
1 - numerator/denominator
}

### this version of nlp_fun did not work on AD; scroll down for updated function

nlp_fun <- function(bin.i, lib_type) {

# bin.i=unique(snp.dt$bin)[10]

### filter
    seqResetFilter(genofile)
    setkey(snp.dt, col)
    seqSetFilter(genofile, variant.id=snp.dt[J(c("missense_variant", "synonymous_variant"))][use==T][bin==bin.i]$id,
                sample.id=samps$sample)

if(lib_type=="individual") {

    ### get dosage matrix
    genomat <- as.matrix(seqGetData(genofile, "$dosage_alt"))
    rownames(genomat) <- seqGetData(genofile, "sample.id")
    colnames(genomat) <- seqGetData(genofile, "variant.id")

    geno.dt <- as.data.table(reshape2::melt(genomat))
    setkey(geno.dt, sample)
    setkey(samps, sample)

    geno.dt <- merge(geno.dt, samps)
    geno.dt

    ### summarize
    geno.geno.dt <- geno.dt[,list(nAA=sum(value==2, na.rm=T), nAa=sum(value==1, na.rm=T), naa=sum(value==0, na.rm=T)),
                            list(pop, variant)]

    geno.geno.dt[,af:=(nAA*2+nAa)/(2*nAA + 2*nAa + 2*naa)]
    geno.geno.dt[,nChr:=(2*nAA + 2*nAa + 2*naa)]



        nlp <- geno.geno.dt[,list(nLocales_poly=sum(af>0 & af<1, na.rm=T),
                                FIS=fis(nHet=sum(nAa, na.rm=T),
                                        nTot=sum(nAA, na.rm=T) + sum(nAa, na.rm=T) + sum(naa, na.rm=T)),
                                tot_nAA=sum(nAA, na.rm=T), tot_nAa=sum(nAa, na.rm=T), tot_naa=sum(naa, na.rm=T),
                                .N,
                                global_af=mean(af, na.rm=T),
                                poly_af=mean(af[af>0 & af<1 & !is.na(af)], na.rm=T)),
                            list(variant=variant)]



    nlp <- merge(nlp, snp.dt, by.x="variant", by.y="id")

    ### Fs

    #try(fst_tmp_new <- seqFst_wc84_exact2(gds=genofile,
    #           sample.id=samps$sample,
    #           population=as.factor(samps$pop),
    #           variant.id=nlp$variant))
    #fst_tmp_new <- as.data.table(fst_tmp_new)
    #setnames(fst_tmp_new, "variant.id", "variant")
    #nlp <- merge(nlp[,-"b"], fst_tmp_new, by="variant")

} else if(lib_type=="pooled") {

    ### getting data
    ad <- seqGetData(genofile, "annotation/format/AD"  )
    dp <- seqGetData(genofile, "annotation/format/DP"  )


    dimnames(ad)$sample <- stringi::stri_escape_unicode(seqGetData(genofile, "sample.id"))
    dimnames(ad)$variant <- seqGetData(genofile, "variant.id")
    dimnames(dp)$sample <- stringi::stri_escape_unicode(seqGetData(genofile, "sample.id"))
    dimnames(dp)$variant <- seqGetData(genofile, "variant.id")

    tmp.ad <-   as.data.table(reshape2::melt(ad))
    tmp.dp <-   as.data.table(reshape2::melt(dp))

    setkey(tmp.ad, sample, variant)
    setkey(tmp.dp, sample, variant)

    m <- merge(tmp.ad, tmp.dp)
    setnames(m, c("value.x", "value.y"), c("ad", "dp"))
    m[,freq:=ad/dp]
    m <- merge(m, samps, by="sample")
    m <- merge(m, snp.dt, by.x="variant", by.y="id")

    ### get FST

    #### format conversion for poolfstat
    #  m[is.na(ad), ad:=0]
    #  m[is.na(dp), dp:=0]

    #  tmp.data.ref.mat <- dcast(data=m[,c("sample", "variant", "ad")], variant~sample, value.var="ad")
    #  tmp.data.dp.mat <- dcast(data=m[,c("sample", "variant", "dp")], variant~sample, value.var="dp")

    #  focalSNPs <- m[,list(.N, freq=mean(freq, na.rm=T)), list(chr, pos, ref, alt)]
    #  setnames(focalSNPs, c("chr", "pos", "ref", "alt"), c("Chromosome", "Position", "RefAllele", "AltAllele"))

    #  tmp.data.summary <- samps[,list(nInd=mean(nInd)), list(sample)]

    #  tmp.pool <- new("pooldata",
    #              npools=dim(tmp.data.summary)[1],
    #              nsnp=dim(tmp.data.dp.mat)[1],
    #              refallele.readcount=as.matrix(tmp.data.ref.mat[,-"variant"]),
    #              readcoverage=as.matrix(tmp.data.dp.mat[,-"variant"]),
    #              poolsizes= tmp.data.summary$nInd * 2,
    #              poolnames = tmp.data.summary$sample,
    #              snp.info = focalSNPs[,c("Chromosome", "Position", "RefAllele", "AltAllele")])

    ## ### calculte pairwise Fst
    #  fst.out <- compute.pairwiseFST(x=tmp.pool, output.snp.values=T, min.maf=0)

        #fst.out <- compute.fstats(x=tmp.pool, computeF3 = F, computeF4 = F, output.pairwise.fst = T)
    #  ### format Fst output
    #    fst.snp.mat <- fst.out@PairwiseSnpFST
    #    rownames(fst.snp.mat) <- focalSNPs$id

    #    fst.snp.long <- as.data.table(reshape2::melt(fst.snp.mat))
    #    fst.snp.long[,samp1:=tstrsplit(Var2, ";")[[1]]]
    #    fst.snp.long[,samp2:=tstrsplit(Var2, ";")[[2]]]
    #    fst.snp.long <- merge(fst.snp.long, focalSNPs, by.x="Var1", by.y="id")
    #    setnames(fst.snp.long, "value", "Fst")

    #  ### format Htot output
    #    Q1.snp.mat <- fst.out@PairwiseSnpQ1
    #    rownames(Q1.snp.mat) <- focalSNPs$id
    #    Q1.snp.long <- as.data.table(reshape2::melt(Q1.snp.mat))

    #  ### format Htot output
    #    Q2.snp.mat <- fst.out@PairwiseSnpQ2
    #    rownames(Q2.snp.mat) <- focalSNPs$id
    #    Q2.snp.long <- as.data.table(reshape2::melt(Q2.snp.mat))

    #  ### combine
    #    qs <- merge(Q1.snp.long, Q2.snp.long, by=c("Var1", "Var2"))
    #    fst.snp.long <- merge(fst.snp.long, qs , by=c("Var1", "Var2"))

    #    fst.snp.long <- fst.snp.long[,-c("Var1", "Var2", "af")]
    #    fst.snp.long[,Htot:=1-value.y]
    #    fst.snp.long[,Hwith:=1-value.x]


    ### summarize
    m.ag <- m[,list(nLocales_poly=length(unique(pop[freq>0 & freq<1 & !is.na(freq)])),
                    global_af=sum(ad, na.rm=T)/sum(dp, na.rm=T),
                    poly_af=mean(freq[freq>0 & freq<1 & !is.na(freq)], na.rm=T)),
                list(variant)]

    m.ag[poly_af>.5, poly_maf:=1-poly_af]
    m.ag[poly_af<.5, poly_maf:=poly_af]

    nlp <- merge(m.ag, snp.dt, by.x="variant", by.y="id")
}

return(nlp)

}


# fixed nlp_fun function for ad$length and ad$data instead of matrix already


mat_to_long <- function(mat) {
  data.table(sample  = rep(rownames(mat), times = ncol(mat)),
             variant = rep(colnames(mat), each  = nrow(mat)),
             value   = as.vector(mat))
}


nlp_fun <- function(bin.i, lib_type) {

# bin.i=unique(snp.dt$bin)[10]

### filter
    seqResetFilter(genofile)
    setkey(snp.dt, col)
    seqSetFilter(genofile, variant.id=snp.dt[J(c("missense_variant", "synonymous_variant"))][use==T][bin==bin.i]$id,
                sample.id=samps$sample)

if(lib_type=="individual") {

    ### get dosage matrix
    genomat <- as.matrix(seqGetData(genofile, "$dosage_alt"))
    rownames(genomat) <- seqGetData(genofile, "sample.id")
    colnames(genomat) <- seqGetData(genofile, "variant.id")

    geno.dt <- as.data.table(reshape2::melt(genomat))
    setkey(geno.dt, sample)
    setkey(samps, sample)

    geno.dt <- merge(geno.dt, samps)
    geno.dt

    ### summarize
    geno.geno.dt <- geno.dt[,list(nAA=sum(value==2, na.rm=T), nAa=sum(value==1, na.rm=T), naa=sum(value==0, na.rm=T)),
                            list(pop, variant)]

    geno.geno.dt[,af:=(nAA*2+nAa)/(2*nAA + 2*nAa + 2*naa)]
    geno.geno.dt[,nChr:=(2*nAA + 2*nAa + 2*naa)]

        nlp <- geno.geno.dt[,list(nLocales_poly=sum(af>0 & af<1, na.rm=T),
                                FIS=fis(nHet=sum(nAa, na.rm=T),
                                        nTot=sum(nAA, na.rm=T) + sum(nAa, na.rm=T) + sum(naa, na.rm=T)),
                                tot_nAA=sum(nAA, na.rm=T), tot_nAa=sum(nAa, na.rm=T), tot_naa=sum(naa, na.rm=T),
                                .N,
                                global_af=mean(af, na.rm=T),
                                poly_af=mean(af[af>0 & af<1 & !is.na(af)], na.rm=T)),
                            list(variant=variant)]

    nlp <- merge(nlp, snp.dt, by.x="variant", by.y="id")

} else if(lib_type=="pooled") {

    ### getting data
    ad <- seqGetData(genofile, "annotation/format/AD")
    dp <- seqGetData(genofile, "annotation/format/DP")

    sample_ids  <- stringi::stri_escape_unicode(seqGetData(genofile, "sample.id"))
    variant_ids <- seqGetData(genofile, "variant.id")

    cum_len  <- cumsum(ad$length)
    has_data <- ad$length == 1
    col_for_variant <- cum_len[has_data]

    ad_mat <- matrix(NA_integer_, nrow=nrow(ad$data), ncol=length(ad$length))
    ad_mat[, has_data] <- ad$data[, col_for_variant]

    dimnames(ad_mat) <- list(sample_ids, variant_ids)
    dimnames(dp)     <- list(sample_ids, variant_ids)

    tmp.ad <- mat_to_long(ad_mat)
    tmp.dp <- mat_to_long(dp)

    tmp.ad[, variant := as.integer(variant)]
    tmp.dp[, variant := as.integer(variant)]

    setkey(tmp.ad, sample, variant)
    setkey(tmp.dp, sample, variant)

    m <- merge(tmp.ad, tmp.dp)
    setnames(m, c("value.x", "value.y"), c("ad", "dp"))
    m[,freq:=ad/dp]
    m <- merge(m, samps, by="sample")
    m <- merge(m, snp.dt, by.x="variant", by.y="id")

    ### summarize
    m.ag <- m[,list(nLocales_poly=length(unique(pop[freq>0 & freq<1 & !is.na(freq)])),
                    global_af=sum(ad, na.rm=T)/sum(dp, na.rm=T),
                    poly_af=mean(freq[freq>0 & freq<1 & !is.na(freq)], na.rm=T)),
                list(variant)]

    m.ag[poly_af>.5, poly_maf:=1-poly_af]
    m.ag[poly_af<.5, poly_maf:=poly_af]

    nlp <- merge(m.ag, snp.dt, by.x="variant", by.y="id")
}

return(nlp)

}