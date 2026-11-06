library(SeqArray)
library(data.table)

mel_gds <- seqOpen("/scratch/ejy4bu/drosophila/DEST_remake/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.snpGap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.gds")
genofile <- mel_gds

snp.dt <- readRDS("/scratch/ejy4bu/drosophila/DEST_remake/snpDT/mel/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.snp_dt_allEffects.rds")

snp.dt[, bin:=floor(variant.id/1000000)]
snp.dt[, col:=effect]
snp.dt[, id := variant.id]
snp.dt[, use := TRUE]


seqResetFilter(genofile)
### test subset

sample_ids <- seqGetData(genofile, "sample.id")
samps <- data.table(
    sample = sample_ids,
    pop = sub("_\\d{4}-\\d{2}-\\d{2}$", "", sample_ids)
)

### Alan's functions
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
    #   seqResetFilter(genofile)
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

        # tmp.ad <-   as.data.table(reshape2::melt(ad))
        # tmp.dp <-   as.data.table(reshape2::melt(dp))

        # setkey(tmp.ad, sample, variant)
        # setkey(tmp.dp, sample, variant)

        # m <- merge(tmp.ad, tmp.dp)
        # setnames(m, c("value.x", "value.y"), c("ad", "dp"))
        # m[,freq:=ad/dp]
        # m <- merge(m, samps, by="sample")

        
        freq <- ad/dp
        m <- as.data.table(reshape2::melt(freq))
        setnames(m, "value", "freq")
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

snp.dt.full <- copy(snp.dt)

### identify bins
bin_ids <- unique(
    snp.dt.full[
        col %in% c("missense_variant", "synonymous_variant") & use == TRUE,
        bin
    ]
)

batch_size <- 100000
all_results <- list()
result_index <- 1

for (bin.i in bin_ids) {

    variant_ids <- unique(
        snp.dt.full[
            col %in% c("missense_variant", "synonymous_variant") &
            use == TRUE & bin == bin.i,
            id
        ]
    )

    batches <- split(
        variant_ids,
        ceiling(seq_along(variant_ids) / batch_size)
    )

    for (i in seq_along(batches)) {

        cat("Bin:", bin.i, "/", length(bin_ids),
            "| Batch:", i, "/", length(batches),
            "| Variants:", length(batches[[i]]), "\n")

        snp.dt <- snp.dt.full[id %in% batches[[i]]]

        all_results[[result_index]] <- nlp_fun(bin.i, "pooled")

        result_index <- result_index + 1

        gc()
    }
    saveRDS(rbindlist(all_results, fill = TRUE), "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/nlp_mel_partial.optimized.rds")
    cat("Saved results through bin:", bin.i, "\n")
}

### restore full snp.dt
snp.dt <- snp.dt.full

### combine results
nlp <- rbindlist(all_results, fill = TRUE)

nrow(nlp)
uniqueN(nlp$variant)

saveRDS(all_results, "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/mel_nlp_complete.optimized.rds")




### test first bin
bin.i <- bin_ids[1]

### get eligible variant IDs
variant_ids <- unique(
    snp.dt.full[
        col %in% c("missense_variant", "synonymous_variant") &
        use == TRUE & bin == bin.i,
        id
    ]
)

### split into batches
batch_size <- 10000
batches <- split(
    variant_ids,
    ceiling(seq_along(variant_ids) / batch_size)
)
length(batches)
### calculate NLP
results <- vector("list", length(batches))

for (i in seq_along(batches)) {
    cat("Batch:", i, "Variants:", length(batches[[i]]), "\n")

    snp.dt <- snp.dt.full[id %in% batches[[i]]]

    results[[i]] <- nlp_fun(bin.i, "pooled")

    gc()
}

### combine results
snp.dt <- snp.dt.full

test_nlp <- rbindlist(results, fill = TRUE)

### end test on first bin