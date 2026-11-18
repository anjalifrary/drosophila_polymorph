### final script to make polyAF

library(SeqArray)
library(data.table)

species <- "mel"    # mel or sim

if (species == "mel") {
    gds_file <- "/scratch/ejy4bu/drosophila/DEST_remake/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.snpGap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.gds"
    snp_file <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/mel/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.snp_dt_allEffects.rds"
    temp_file <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/nlp/nlp_mel_partial.16_11_2026.rds"
    nlp_file <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/nlp/mel_nlp.16_11_2026.rds"
} else if (species == "sim"){
    gds_file <- "/scratch/ejy4bu/drosophila/DEST_remake/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.norm.snpGap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.merged.biallelic.dm6.sorted.m2.M2.gds"
    snp_file <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/sim/dest.sim.all.SNAPE.001.50.20Nov2025_sim.snp_dt_allEffects.10_2026.rds"
    temp_file <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/nlp/nlp_sim_partial.16_11_2026.rds"
    nlp_file <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/nlp/sim_nlp.16_11_2026.rds"
} else {
    stop("invalid species")
}

genofile <- seqOpen(gds_file)
snp.dt <- readRDS(snp_file)
seqResetFilter(genofile)



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
    if(whichMean=="polymorphic") return(mean(tmp.dt.ag[freq>0 & freq<1]$freq, na.rm=TRUE))
    if(whichMean=="first") return(tmp.dt.ag$freq[1])
    if(whichMean=="US_Nor_Ral") return(tmp.dt.ag[pop=="US_Nor_Ral"]$freq)

  }

nlp_fun <- function(variant_ids, snp_dt_batch, lib_type = "pooled") {

    # bin.i=unique(snp.dt$bin)[10]

    ### which SNP sets
    # setkey(snp.dt, col)

    # if(colTypes=="NS_SYN") {
    #     ids2use <- snp.dt[J(c("missense_variant", "synonymous_variant"))][use==T][bin==bin.i]$id
    # } else if(colTypes=="all") {
    #     ids2use <- snp.dt[use==T][bin==bin.i]$id
    # }

    ### filter
    seqResetFilter(genofile)
    #   setkey(snp.dt, col)
    #   variant_ids <- snp.dt[J(c("missense_variant", "synonymous_variant"))][use==TRUE & bin == bin.i, unique(id)]
    seqSetFilter(genofile, variant.id=variant_ids, sample.id=samps$sample)
    # cat("Variants selected:", length(seqGetData(genofile, "variant.id")), "\n")
    
    selected_ids <- seqGetData(genofile, "variant.id")
    cat("Variants selected:", length(selected_ids), "\n")

    stopifnot(
        length(selected_ids) == length(variant_ids),
        setequal(selected_ids, variant_ids)
    )

    if(lib_type=="individual") {

      ### get dosage matrix
        genomat <- as.matrix(seqGetData(genofile, "$dosage_alt"))
        rownames(genomat) <- seqGetData(genofile, "sample.id")
        colnames(genomat) <- seqGetData(genofile, "variant.id")

        geno.dt <- as.data.table(reshape2::melt(genomat))
        # setkey(geno.dt, sample)
        # setkey(samps, sample)

        geno.dt <- merge(geno.dt, samps, by = "sample")
        # geno.dt

      ### summarize
        geno.geno.dt <- geno.dt[, .(nAA=sum(value==2, na.rm=T), 
                                        nAa=sum(value==1, na.rm=T), 
                                        naa=sum(value==0, na.rm=T)),
                                by = .(pop, variant)]

        geno.geno.dt[,af:=(2*nAA+nAa)/(2*nAA + 2*nAa + 2*naa)]
        geno.geno.dt[,nChr:=(2*nAA + 2*nAa + 2*naa)]

        nlp <- geno.geno.dt[, .(nLocales_poly=sum(af>0 & af<1, na.rm=T),
                                   nLocales_total=sum(!is.na(af)), 
                                   nLocales_missing=sum(is.na(af)),
                                   FIS=fis(nHet=sum(nAa, na.rm=T),
                                           nTot=sum(nAA, na.rm=T) + sum(nAa, na.rm=T) + sum(naa, na.rm=T)),
                                   tot_nAA=sum(nAA, na.rm=T), tot_nAa=sum(nAa, na.rm=T), tot_naa=sum(naa, na.rm=T),
                                   global_af=mean_allele_freq_fun(freq=af, pop=pop, "global"),                   
                                   poly_af=mean_allele_freq_fun(freq=af, pop=pop, "polymorphic")
                                #    , randomPop_af=mean_allele_freq_fun(freq=af, pop=pop, "first")
                                ), by = .(variant)]



        nlp <- merge(nlp, snp_dt_batch, by.x="variant", by.y="variant.id")

    } else if(lib_type=="pooled") {

      ### getting data
      if (length(seqGetData(genofile, "variant.id")) == 0) {
            stop("Zero variants selected")
        }
        ad <- seqGetData(genofile, "annotation/format/AD"  )
        dp <- seqGetData(genofile, "annotation/format/DP"  )


        if (class(ad)[1] == "SeqVarDataList") {
            ### sim only:
            # Confirm that each variant has either zero or one AD value
            if (!all(ad$length %in% c(0, 1))) {
                stop("Unexpected AD structure: some variants have multiple AD values.")
            } 

            ad_variant <- rep(seq_along(ad$length), ad$length)
            ad_full <- matrix(NA_integer_, nrow = nrow(dp), ncol = length(ad$length))
            ad_full[, ad_variant] <- ad$data
        } else {
            ### mel only
            ad_full <- ad
        }
        dimnames(ad_full) <- list(
            sample = stringi::stri_escape_unicode(seqGetData(genofile, "sample.id")),
            variant = seqGetData(genofile, "variant.id")
        )
        dimnames(dp) <- list(
            sample = stringi::stri_escape_unicode(seqGetData(genofile, "sample.id")),
            variant = seqGetData(genofile, "variant.id")
        )
        tmp.ad <- as.data.table(reshape2::melt(ad_full))
        tmp.dp <- as.data.table(reshape2::melt(dp))        

        # # ### mel only:
        # dimnames(ad)$sample <- stringi::stri_escape_unicode(seqGetData(genofile, "sample.id"))
        # dimnames(ad)$variant <- seqGetData(genofile, "variant.id")

        # dimnames(dp)$sample <- stringi::stri_escape_unicode(seqGetData(genofile, "sample.id"))
        # dimnames(dp)$variant <- seqGetData(genofile, "variant.id")

        # tmp.ad <-   as.data.table(reshape2::melt(ad))
        # tmp.dp <-   as.data.table(reshape2::melt(dp))
        # # ### end mel only ^

        tmp.ad[, variant := as.integer(as.character(variant))]
        tmp.dp[, variant := as.integer(as.character(variant))]
        
        setkey(tmp.ad, sample, variant)
        setkey(tmp.dp, sample, variant)

        m <- merge(tmp.ad, tmp.dp, by = c("sample", "variant"))
        setnames(m, c("value.x", "value.y"), c("ad", "dp"))
        m[, freq:=ad/dp]
        m[is.na(ad) | is.na(dp) | dp <= 0, freq := NA_real_]
        m <- merge(m, samps, by="sample")
        # m <- merge(m, snp.dt, by.x="variant", by.y="variant.id")
  
        ### collapse by locality
        # m2 <- m[,list(freq=mean(freq, na.rm=T)), list(variant, pop)]
        
        
      ### summarize
        # if(spp%like%"Drosophila_melanogaster") {
        #     randomPop_pop="US_Nor_Ral"
        # } else {
        #     randomPop_pop="first"
        # }

        m.ag <- m[, .(nLocales_poly=length(unique(pop[freq>0 & freq<1 & !is.na(freq)])),
                        nLocales_total=length(unique(pop[!is.na(freq)])), 
                        nLocales_missing=length(unique(pop[is.na(freq)])),
                        global_af=mean_allele_freq_fun(freq=freq, pop=pop, "global"),                    ### hold over from an older version: sum(ad, na.rm=T)/sum(dp, na.rm=T),
                        poly_af=mean_allele_freq_fun(freq=freq, pop=pop, "polymorphic")
                        # ,
                        # randomPop_af=mean_allele_freq_fun(freq=freq, pop=pop, "first")
                        ), by = .(variant)]

        nlp <- merge(m.ag, snp_dt_batch, by.x = "variant", by.y = "variant.id")

        # nlp <- merge(m.ag, snp.dt, by.x="variant", by.y="variant.id")
        
    }
    stopifnot(
        !anyDuplicated(nlp$variant),
        all(nlp$variant %in% variant_ids)
    )

    return(nlp)

}

snp.dt[, bin:=floor(variant.id/1000000)]
snp.dt[, col:=effect]
snp.dt[, use := TRUE]

sample_ids <- seqGetData(genofile, "sample.id")
samps <- data.table(
    sample = sample_ids,
    pop = sub("_\\d{4}-\\d{2}-\\d{2}$", "", sample_ids)
)

snp.dt.full <- copy(snp.dt)

### identify bins
bin_ids <- unique(
    snp.dt.full[
        col %in% c("missense_variant", "synonymous_variant") & use == TRUE,
        bin
    ]
)

batch_size <- 100000
nlp_results <- data.table()
nlp_results <- readRDS("/scratch/ejy4bu/drosophila/DEST_remake/snpDT/nlp/nlp_mel_partial.16_11_2026.rds")

for (bin.i in bin_ids) {

    variant_ids <- unique(
        snp.dt.full[
            col %in% c("missense_variant", "synonymous_variant") &
            use == TRUE & bin == bin.i,
            variant.id])

    batches <- split(variant_ids, ceiling(seq_along(variant_ids) / batch_size)
    )

    for (i in seq_along(batches)) {
        batch_ids <- batches[[i]]
        batch_ids <- setdiff(batch_ids, nlp_results$variant)
        if (length(batch_ids) == 0) next

        cat("Bin:", bin.i, "/", length(bin_ids),
            "| Batch:", i, "/", length(batches),
            "| Variants:", length(batch_ids), "\n")

        snp.dt.batch <- snp.dt.full[variant.id %in% batch_ids]
        stopifnot(
            nrow(snp.dt.batch) == length(batch_ids),
            uniqueN(snp.dt.batch$variant.id) == length(batch_ids),
            setequal(snp.dt.batch$variant.id, batch_ids)
        )

        # batch_result <- nlp_fun(batches[[i]], "pooled")
        batch_result <- nlp_fun( variant_ids = batch_ids, snp_dt_batch = snp.dt.batch, lib_type = "pooled")

        cat("Batch variants:", length(batch_ids), "\n")
        cat("SNP table rows:", nrow(snp.dt.batch), "\n")
        cat("NLP result rows:", nrow(batch_result), "\n")
        
        stopifnot(
            nrow(batch_result) == length(batch_ids),
            uniqueN(batch_result$variant) == length(batch_ids),
            setequal(batch_result$variant, batch_ids)
        )

        nlp_results <- rbindlist(
            list(nlp_results, batch_result),
            fill = TRUE
        )
         stopifnot( !anyDuplicated(nlp_results$variant) )

        saveRDS(nlp_results, temp_file)

        cat("Saved", nrow(nlp_results), " variants through bin:", bin.i, ", batch: ", i, "\n")
        rm(batch_result, snp.dt.batch, batch_ids)
        gc()
    }
    rm(variant_ids, batches)
    gc()
}

### restore full snp.dt
snp.dt <- snp.dt.full
stopifnot(
    !anyDuplicated(nlp_results$variant),
    uniqueN(nlp_results$variant) == uniqueN(
        snp.dt.full[
            col %in% c("missense_variant", "synonymous_variant") & use == TRUE,
            variant.id])
)

### combine results; filter to only variants polymorphic in at least 1 population
nlp <- nlp_results
nlp <- nlp[nLocales_poly > 0]

nlp.miss <- nlp[, .(.N), by = nLocales_total]
nlp.miss <- nlp.miss[order(-nLocales_total)]
nlp.miss[, csN := cumsum(N)]
nlp.miss[, fcsN := csN / max(nlp.miss$csN)]

projTarget <- nlp.miss[which.min(abs(fcsN - 0.95)), nLocales_total]
cat("Projection target:", projTarget, "populations\n")

nlp[
    nLocales_total >= projTarget,
    nLocales_poly_proj_round :=
    round((nLocales_poly / nLocales_total) * projTarget)
]

nlp[
    nLocales_total >= projTarget,
    nLocales_poly_proj_binomial :=
    rbinom(
        # dim(nlp[nLocales_total >= projTarget])[1],
        .N,
        size = projTarget,
        prob = (nLocales_poly / nLocales_total)
    )
]
stopifnot(
    !anyDuplicated(nlp$variant),
    all(nlp$nLocales_poly > 0)
)

nrow(nlp)
uniqueN(nlp$variant)
cat("Final polymorphic variants:", nrow(nlp), "\n")
cat("Unique variants:", uniqueN(nlp$variant), "\n")
table(is.na(nlp$nLocales_poly_proj_round), nlp$nLocales_total)


# saveRDS(all_results, "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/mel_nlp_complete.optimized.rds")
# saveRDS(nlp, "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/mel_nlp_complete.new.NArm.rds")
saveRDS(nlp, nlp_file)
cat("Saved to file: ", nlp_file, "\n")



### test
snp.dt.full <- copy(snp.dt)

# Select the first 10 variants that will actually be analyzed
test_ids <- unique(
    snp.dt.full[
        col %in% c("missense_variant", "synonymous_variant") &
        use == TRUE,
        id
    ]
)[1:10]

cat("Test IDs:\n")
print(test_ids)

cat("\nID class in snp.dt:\n")
print(class(snp.dt.full$id))

cat("\nNumber of unique test IDs:\n")
print(uniqueN(test_ids))


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
    cat("Bin:", bin.i, "/", length(bin_ids),
            "| Batch:", i, "/", length(batches),
            "| Variants:", length(batches[[i]]), "\n")
    snp.dt <- snp.dt.full[id %in% batches[[i]]]

    results[[i]] <- nlp_fun(bin.i, "pooled")

    gc()
}

### combine results
snp.dt <- snp.dt.full

test_nlp <- rbindlist(results, fill = TRUE)

### end test on first bin