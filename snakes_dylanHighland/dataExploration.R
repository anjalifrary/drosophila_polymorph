library(data.table)
library(SeqArray)

tigris <- "/scratch/ejy4bu/dylan_snakes/snakes/data/gds/crotalus_tigris.snps.NS_SYN.annotate.gds"
steph <- "/scratch/ejy4bu/dylan_snakes/snakes/data/gds/crotalus_stephensi.snps.NS_SYN.annotate.gds"
pyrrhus <- "/scratch/ejy4bu/dylan_snakes/snakes/data/gds/crotalus_pyrrhus.snps.NS_SYN.annotate.gds"
mitch <- "/scratch/ejy4bu/dylan_snakes/snakes/data/gds/crotalus_mitchellii.snps.NS_SYN.annotate.gds"

snp.dt <- foreach(gds.fn=fn)%do%{
    gds.fn <- tigris
    genofile <- seqOpen(gds.fn)
    seqResetFilter(genofile)

    dt <- data.table(
        chr = seqGetData(genofile, "chromosome"),
        pos = seqGetData(genofile, "position"), 
        nAlleles = seqGetData(genofile, "$num_allele"), 
        id = seqGetData(genofile, "variant.id"),
        af = seqGetData(genofile, "annotation/info/AF"), 
        alleles = seqGetData(genofile, "allele"),
        ref = sub(",.*", "", alleles),
        alt = sub(".*,", "", alleles),
        ref2=seqGetData(genofile, "$ref"),
        alt2=seqGetData(genofile, "$alt")    
    )

    dt[, is_snp := (nchar(ref) == 1 & nchar(alt) == 1)]
    dt[, is_indel := (nchar(ref) != nchar(alt))]
}


is_snp <- nchar(ref) == 1 & nchar(alt) == 1
    is_indel <- nchar(ref) != nchar(alt)

isSNP <- function(ref, alt) {
    ## TRUE if ref and every alt allele are single nucleotides (A/C/G/T).
    ## Handles multiallelic alts ("A,T"); indels, MNPs, "*", "<DEL>", "." etc. return FALSE.

        ref <- toupper(ref)
        alt <- toupper(alt)
        refOK <- grepl("^[ACGT]$", ref)
        altOK <- grepl("^[ACGT](,[ACGT])*$", alt)
        refOK & altOK
    }


dt <- data.table(
    chr = seqGetData(gds_file, "chromosome"), 
    pos = seqGetData(gds_file, "position"),
    nAlleles = seqGetData(gds_file, "$num_allele"),
    id = seqGetData(gds_file, "variant.id"),
    af = seqGetData(gds_file, "annotation/info/AF") 
    # ,
    # maf = seqGetData(gds_file, "annotation/info/MAF"),
    # n_samps = seqGetData(gds_file, "annotation/info/NS")
)


seqClose(tigris)
seqClose(steph)
seqClose(pyrrhus)
seqClose(mitch)