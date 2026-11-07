### libraries
    library(data.table)
    library(SeqArray)
    library(foreach)

### convert
    fn <- list.files("/Users/alanbergland/Documents/GitHub/misc/snakes/data", "vcf.gz", full.name=T)

    foreach(vcf.fn=fn)%do%{
        # vcf.fn <- fn[1]
        gds.fn=gsub(".vcf.gz", ".gds", vcf.fn)

        seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", verbose=T, optimize=T)
    }

