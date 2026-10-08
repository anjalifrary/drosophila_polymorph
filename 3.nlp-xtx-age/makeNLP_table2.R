### FROM ALAN 10/8/26

# look at lines 18/19 for poly & global AF

### common functions

  fis <- function(nHet, nTot, p=.5) {
    numerator <- (nHet/nTot)
    denominator <- (nTot/(nTot-1))*(2*p*(1-p) - (nHet/nTot)/(2*nTot))
    1 - numerator/denominator
  }
  
  mean_allele_freq_fun <- function(freq, pop, whichMean) {
    # freq <- m[variant==100001]$freq; pop <- m[variant==100001]$pop
    tmp.dt <- data.table(freq=freq, pop=pop)
    tmp.dt.ag <- tmp.dt[,list(freq=mean(freq, na.rm=T)), list(pop)]

    if(whichMean=="global") return(mean(tmp.dt.ag$freq, na.rm=T))
    if(whichMean=="polymorphic") return(mean(tmp.dt.ag[freq>0 & freq<1]$freq))
    if(whichMean=="first") return(tmp.dt.ag$freq[1])
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

        geno.geno.dt[,af:=(2*nAA+nAa)/(2*nAA + 2*nAa + 2*naa)]
        geno.geno.dt[,nChr:=(2*nAA + 2*nAa + 2*naa)]

         nlp <- geno.geno.dt[,list(nLocales_poly=sum(af>0 & af<1, na.rm=T),
                                   nLocales_total=sum(!is.na(af)), 
                                   nLocales_missing=sum(is.na(af)),
                                   FIS=fis(nHet=sum(nAa, na.rm=T),
                                           nTot=sum(nAA, na.rm=T) + sum(nAa, na.rm=T) + sum(naa, na.rm=T)),
                                   tot_nAA=sum(nAA, na.rm=T), tot_nAa=sum(nAa, na.rm=T), tot_naa=sum(naa, na.rm=T),
                                   global_af=mean_allele_freq_fun(freq=af, pop=pop, "global"),                   
                                   poly_af=mean_allele_freq_fun(freq=af, pop=pop, "polymorphic"),
                                   randomPop_af=mean_allele_freq_fun(freq=af, pop=pop, "first")),
                              list(variant=variant)]



        nlp <- merge(nlp, snp.dt, by.x="variant", by.y="id")

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
  
        
      ### summarize

        m.ag <- m[,list(nLocales_poly=length(unique(pop[freq>0 & freq<1 & !is.na(freq)])),
                        nLocales_total=length(unique(pop[!is.na(freq)])), 
                        nLocales_missing=length(unique(pop[is.na(freq)])),
                        global_af=mean_allele_freq_fun(freq=freq, pop=pop, "global"),                    ### hold over from an older version: sum(ad, na.rm=T)/sum(dp, na.rm=T),
                        poly_af=mean_allele_freq_fun(freq=freq, pop=pop, "polymorphic"),
                        randomPop_af=mean_allele_freq_fun(freq=freq, pop=pop, "first")),
                  list(variant)]

        nlp <- merge(m.ag, snp.dt, by.x="variant", by.y="id")
    }

    return(nlp)

  }