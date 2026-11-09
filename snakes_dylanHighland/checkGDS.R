### libraries
    library(data.table)
    library(SeqArray)
    library(foreach)
    options(width = 250)

### functions
    # source("/Users/alanbergland/Documents/GitHub/misc/snakes/scripts/functions.R")

### convert
    # fn <- list.files("/Users/alanbergland/Documents/GitHub/misc/snakes/data", "gds", full.name=T)

### make snp.dts
    snp.dt <- foreach(gds.fn=fn)%do%{
        # gds.fn <- fn[1]
        genofile <- seqOpen(gds.fn)
        seqResetFilter(genofile)
      
        snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       nAlleles=seqGetData(genofile, "$num_allele"),
                       id=seqGetData(genofile, "variant.id"),
                       ref=seqGetData(genofile, "$ref"),
                       alt=seqGetData(genofile, "$alt"))
        snp.dt[,isSNP:=isSNP(ref, alt)]
        
        snp.dt[,spp:=gsub(".snps.NS", "", tstrsplit(gsub(".*data/", "", gds.fn), "_")[[2]])]
        snp.dt[,gdsFN:=gsub(".*data/", "", gds.fn)]
        
        snp.dt[,ref_0_freq:=seqAlleleFreq(genofile)]
        snp.dt[,alt_1_freq:=seqAlleleFreq(genofile, ref.allele=1L)]
        snp.dt[,alt_2_freq:=seqAlleleFreq(genofile, ref.allele=2L)]
        
        snp.dt[, c("ref_used", "alt_used", "ref_freq", "alt_freq") :=
               getBiallelic(ref, alt, ref_0_freq, alt_1_freq, alt_2_freq)]

        eff.dt <- getEFF(genofile)
        setorder(eff.dt, id, eff_allele, ref_codon, na.last=TRUE)
        eff.dt <- eff.dt[, .SD[1], .(id, eff_allele)]
        snp.dt <- merge(snp.dt, eff.dt, by.x=c("id", "alt_used"), by.y=c("id", "eff_allele"),
                        all.x=TRUE, sort=FALSE)


        seqClose(genofile)
        return(snp.dt)
    }   

    snp.dt <- rbindlist(snp.dt)
    snp.dt <- snp.dt[isSNP==T]

### check dimensions
    snp.dt[,list(.N), list(spp)]

### get number of polymorphic species at the SNP level
    snp.dt[,spp_poly:=sum(ref_freq > 0 & ref_freq < 1, na.rm=T), list(chr, pos)]
    snp.dt[!is.na(aa_pos), codon_poly:=sum(ref_freq > 0 & ref_freq < 1, na.rm=T), list(spp,transcript, aa_pos)][!is.na(aa_pos), codon_poly_max:=max(codon_poly), list(transcript, aa_pos)]

### look at example of TSP NS
    snp.dt[pos==17238384][ref_freq>0 & ref_freq<1]

### convergent?
    poly.pairs <- classifyPoly(snp.dt)
    poly.pairs[, .N, list(class, type)]

### look at example of convergent NS
    poly.pairs[class=="convergent"]
    snp.dt[pos==17464684][ref_freq>0 & ref_freq<1]

### look at example of "different"
    poly.pairs[class=="different"]
    snp.dt[aa_pos==135 & gene=="venom_SVMP-1-mdc"][pos==17336109][,c("chr", "pos", "spp", "ref_freq", "alt_freq", "ref_codon", "alt_codon", "ref_aa", "alt_aa", "codon_poly")]
    snp.dt[aa_pos==135 & gene=="venom_SVMP-1-mdc"][pos==17336111][,c("chr", "pos", "spp", "ref_freq", "alt_freq", "ref_codon", "alt_codon", "ref_aa", "alt_aa", "codon_poly")]
    poly.pairs[class=="different"][aa_pos==306 & gene=="venom_SVMP-9-mad"]

### make wide table
    snp.dt.w <- dcast(snp.dt, chr+pos+effect+gene~spp, value.var=c("ref_codon", "alt_codon", "ref_aa", "alt_aa"))
    snp.dt.w
    