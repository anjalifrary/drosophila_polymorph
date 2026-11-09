library(data.table)
library(SeqArray)

tigris <- seqOpen("/scratch/ejy4bu/dylan_snakes/snakes/data/crotalus_tigris.snps.NS_SYN.annotate.gds")
steph <- seqOpen("/scratch/ejy4bu/dylan_snakes/snakes/data/crotalus_stephensi.snps.NS_SYN.annotate.gds")
pyrrhus <- seqOpen("/scratch/ejy4bu/dylan_snakes/snakes/data/crotalus_pyrrhus.snps.NS_SYN.annotate.gds")
mitch <- seqOpen("/scratch/ejy4bu/dylan_snakes/snakes/data/crotalus_mitchellii.snps.NS_SYN.annotate.gds")


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
