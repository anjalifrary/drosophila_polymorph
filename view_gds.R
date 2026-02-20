# .libPaths(c("~/Rlibs", .libPaths()))
# BiocManager::install("SeqArray" force=T)

library(SeqArray)
# library(SeqVarTools)

library(data.table)
library(ggplot2)
library(patchwork)
library(foreach)
library(doMC)
registerDoMC(10)


out_file <- "/scratch/ejy4bu/drosophila/gds_analysis/view_gds.txt"
if(!file.exists(out_file)) file.create(out_file)

# load melanogaster gds file
mel_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.gds"
mel_gds <- seqOpen(mel_file)
mel_gds

# load simulans gds file
sim_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds"
sim_gds <- seqOpen(sim_file)
sim_gds 

genofile <- mel_gds

# # reference allele freq of each variant
# mel_freq <- seqAlleleFreq(mel_gds)
# head(mel_freq)
# summary(mel_freq)
# mel_df <- data.frame(freq = mel_freq)

# sim_freq <- seqAlleleFreq(sim_gds)
# head(sim_freq)
# summary(sim_freq)
# sim_df <- data.frame(freq = sim_freq)

# # take a subset of variants by id 
# seqSetFilter(mel_gds, variant.id=100:105)
# seqGetData(mel_gds,"annotation/info/ANN")

seqResetFilter(genofile)
snp.dt <- data.table(
    chr=seqGetData(genofile, "chromosome"),
    pos=seqGetData(genofile, "position"),
    nAlleles=seqGetData(genofile, "$num_allele"),
    id=seqGetData(genofile,"variant.id"))

snp.dt <- snp.dt[nAlleles==2] ##subset to sites with 2 alleles
seqSetFilter(genofile, variant.id = snp.dt$id)

# pull global allele freq field from gds
snp.dt[,af:=seqGetData(genofile,"annotation/info/AF")$data] # global average frequency

# load metadata?
in_dir <- "/project/berglandlab/anjali/drosophila_polymorphism/metadata"
mel_meta <- fread(file.path(in_dir,"dest_v2.samps_24Aug2024.xa.csv"))
samps <- mel_meta

### load this function
getData <- function(snps=snp.dt[pos==14617051 & chr=="2L"], samples=samps) {
    # snps=snp.dt[pos==14617051 & chr=="2L"]; samples=samps
    # this position = adh polymorphism? 

    ### filter to target
    seqSetFilter(genofile, variant.id=snps$id)

    ### get annotations
    message("Annotations")
    tmp <- seqGetData(genofile, "annotation/info/ANN")
    len1 <- tmp$length
    len2 <- tmp$data

    snp.dt1 <- data.table(len=rep(len1, times=len1),
                          ann=len2,
                          id=rep(snps$id, times=len1))

    # Extract data between the 2nd and third | symbol
    snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
    snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]

    # Collapse additional annotations to original SNP vector length
    snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), gene=paste(gene, collapse=",")),
                          list(variant.id=id)]

    snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
    snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]

    ### get frequencies
    message("Allele Freqs")

    ad <- seqGetData(genofile, "annotation/format/AD")
    dp <- seqGetData(genofile, "annotation/format/DP")

    # ad_mat <- ad$data
    # dp_mat <- dp$data

    n_samples <- length(seqGetData(genofile, "sample.id"))
    n_variants <- length(seqGetData(genofile, "variant.id"))

    if(class(dp)[1]!="SeqVarDataList") {
      dp <- list(data=dp)
    }
    if(class(ad)[1]!="SeqVarDataList") {
      ad <- list(data=ad)
    }

    af <- data.table(ad=as.vector(ad$data),
                     dp=as.vector(dp$data),
                     sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                     variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))

    # af <- data.table(
    #     ad = as.vector(ad_mat),
    #     dp = as.vector(dp_mat),
    #     sampleId = rep(seqGetData(genofile, "sample.id"), times= n_variants),
    #     variant.id = rep(seqGetData(genofile, "variant.id"), each= n_samples)
    # )

    ### tack them together
    message("merge")
    afi <- merge(af, snp.dt1.an, by="variant.id")
    afi <- merge(afi, snp.dt, by.x="variant.id", by.y="id")
    afi[,af:=ad/dp] # calculate average allele frequency 

    ### calculate effective read-depth
    afis <- merge(afi, samples[,c("sampleId", "nFlies", "locality",
                                  "lat", "long", "continent", "country", "province", "city",
                                  "min_day", "max_day", "min_month", "max_month", "year", "jday",
                                  "bio_rep", "tech_rep", "exp_rep", "loc_rep", "subsample", "sampling_strategy",
                                  "SRA_Accession"), with=F], by="sampleId")

    # effective allele frequency                  
    afis[chr=="X|Y", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
    afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
    afis[,af_nEff:=round(af*nEff)/nEff]

    setnames(afis, "col", "annotation")
    afis[,-c("n"), with=F]
}

# example
result <- getData(
    snps = snp.dt[pos==14617051 & chr=="2L"],
    samples = samps
)
result[, .N, by = annotation]

# result_region



### histogram of allele frequency

### allele freq = (# alt alleles observed)/(# total alleles observed)
# near 1 => fixed alleles
# near 0.5 => SNP
# near 0 => rare allele

library(ggplot2)
# mel:
p <- ggplot(mel_df, aes(x=freq)) + 
    geom_histogram(binwidth=0.05) +
    labs(title = "Distribution of Mel Allele Frequencies",
    x = "Allele Frequency",
    y = "Count")
ggsave("/scratch/ejy4bu/drosophila/gds_analysis/mel_freq_histogram.jpg", plot = p, width=8,height=6)

# sim
p <- ggplot(sim_df, aes(x=freq)) + 
    geom_histogram(binwidth=0.05) +
    labs(title = "Distribution of Sim Allele Frequencies",
    x = "Allele Frequency",
    y = "Count")
ggsave("/scratch/ejy4bu/drosophila/gds_analysis/sim_freq_histogram.jpg", plot = p, width=8,height=6)
