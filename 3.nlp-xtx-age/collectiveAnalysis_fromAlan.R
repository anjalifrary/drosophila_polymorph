# ijob -A berglandlab -c20 -p standard --mem=200G
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

### libraries
library(SeqArray)
library(data.table)
#library(ggplot2)
library(foreach)
library(doMC)
registerDoMC(20)
library(readxl)
library(SNPRelate)
library(poolfstat)

output_directory <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/"

spp <- "Drosophila_simulans"
# spp <- "Drosophila_melanogaster"

# relevant section for anjali

# test <- nlp_fun(bin.i=bin.i, lib_type="pooled")
# str(test)

nlp <- foreach(bin.i=unique(snp.dt$bin), .errorhandling="remove") %dopar% {
  nlp_fun(bin.i=bin.i, lib_type="pooled")
}
nlp <- rbindlist(nlp)
nlp <- nlp[nLocales_poly>0]

# these should match:
length(unique(snp.dt$bin))
length(unique(nlp$bin))

save(nlp, file=paste(output_directory, spp, ".", "17_06_2026", ".nlpTable.Rdata", sep=""))



#############################
# SKIP THIS STUFF FOR NOW
#############################

### which species to work on
args <- commandArgs(trailingOnly = TRUE)
sppNum <- as.numeric(args[1])

#sppNum <- 19

### define directory


# ### generate meta-meta object
# source("~/misc/1000G/done/1.collectiveAnalysis/0.createMetadata.R")
# names(species)[sppNum]
# spp <- species[[sppNum]]$spp

# ### open GDS file
# try(genofile <- seqOpen(species[[sppNum]]$gdsFile, allow.duplicate=TRUE))

# ### generate metadata file
# source(species[[sppNum]]$sample_metadata)
# str(samps)

# ### what species?
# message(paste("working on species: ", spp))

######################
### make SNP table ###
######################
### version
SNP_table_build_date <- "17_06_2026"  ### initial attemp; DMY

# ### run or load
# if(!file.exists(paste(output_directory, spp, ".", SNP_table_build_date, ".snpTable.Rdata", sep=""))) {
#     source("~/misc/1000G/done/1.collectiveAnalysis/1.makeSNP_table.R")
#     snp.dt[,spp:=spp]
#     save(snp.dt, file=paste(output_directory, spp, ".", SNP_table_build_date, ".snpTable.Rdata", sep=""))
# } else {
#     message("Current SNP table already exists; loading")
#     load(paste(output_directory, spp, ".", SNP_table_build_date, ".snpTable.Rdata", sep=""))
# }

### inspect
str(snp.dt)
setkey(snp.dt, col)

######################
### make NLP table ###
######################
### version
NLP_table_build_date <- "10_03_2026"  ### initial attemp; DMY

spp <- Drosophila_simulans
 

### run or load
source("~/misc/1000G/done/1.collectiveAnalysis/2.makeNLP_table.R")

if(!file.exists(paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.Rdata", sep=""))) {

    nlp <- foreach(bin.i=unique(snp.dt$bin), .errorhandling="remove")%dopar%{
    # bin.i <- 10
    nlp_fun(bin.i=bin.i, lib_type=species[[sppNum]]$libraryType)
    }
    nlp <- rbindlist(nlp)
    nlp <- nlp[nLocales_poly>0]

    message(paste("saving NLP dimensions:", dim(nlp), sep=" "))
    save(nlp, file=paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.Rdata", sep=""))

} else {
    message("Current basic NLP table already exists; loading")
    load(file=paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.Rdata", sep=""))
}
nlp.orig <- nlp

### inspect
str(nlp)

###############################
### load in ParaMask output ###
###############################
if(!file.exists(paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.paramask.Rdata", sep=""))) {
message("compiling paramask or other CNVs")

if(!is.na(species[[sppNum]]$paramask)) {
    pm <- fread(species[[sppNum]]$paramask)
    setnames(pm, c("Chromosome", "Position"), c("chr", "pos"))

    if(species[[sppNum]]$spp=="Crotalus_helleri") pm[,chr:=gsub("chr", "", chr)]
    if(species[[sppNum]]$spp=="Mimulus_gutattus") pm[,chr:=as.character(chr)]

    setkey(pm, chr, pos)
    setkey(nlp, chr, pos)

    nlp.pm <- merge(nlp, pm, all.x=T)
    try(nlp <- nlp.pm)

}

#########################################
### other CNV calls; somewhat bespoke ###
#########################################
    if(species[[sppNum]]$spp=="Arabidopsis_thaliana") {
    cnv <- as.data.table(read_excel("~/misc/1000G/done/arabidopsis/TPC2019-LSB-00640R2_Supplemental_Data_Set_1_8.xlsx", skip=4, sheet="Supplemental Data Set 1"))
    cnv[,Chr:=gsub("Chr", "", Chr)]
    setnames(cnv, "Chr", "chr")
    nlp[,Start:=pos]
    nlp[,Stop:=pos+1]

    setkey(cnv, chr, Start, Stop)
    setkey(nlp, chr, Start, Stop)


    nlp.cnv <- foverlaps(nlp, cnv[,c("chr", "Start", "Stop", "CNV_ID")], nomatch=NA)

    #str(nlp.cnv)
    nlp.cnv <- nlp.cnv[,-c("Start", "Stop", "i.Start", "i.Stop")]
    nlp.cnv[,finalClass:=as.numeric(is.na(CNV_ID))]
    #str(nlp.cnv)
    nlp <- nlp.cnv
    }

message(paste("saving NLP + paramask/cnv dimensions:", dim(nlp), sep=" "))
save(nlp, file=paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.paramask.Rdata", sep=""))

} else {
message("loaing paramask nlp")

load(file=paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.paramask.Rdata", sep=""))
}

##########################
### mappability scores ###
##########################
if(!file.exists(paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.paramask.genmap.Rdata", sep=""))) {
message("tacking in genmap")

if(!is.na(species[[sppNum]]$genmap)) {
    genmap <- fread(species[[sppNum]]$genmap)

    setnames(genmap, c("V1", "V2", "V3"), c("chr", "start", "stop"))

    if(species[[sppNum]]$spp=="Homo_sapiens" & !(any(nlp$chr%like%"chr"))){
    nlp[,chr:=paste("chr", chr, sep="")]
    }
    if(species[[sppNum]]$spp%like%"Anopheles") {
    genmap[chr=="NC_004818.2", chr:= "X" ]
    genmap[chr=="NT_078265.2", chr:= "2L"]
    genmap[chr=="NT_078266.2", chr:= "2R"]
    genmap[chr=="NT_078267.5", chr:= "3L"]
    genmap[chr=="NT_078268.4", chr:= "3R"]
    }
    if(species[[sppNum]]$spp%like%"Arabidopsis_thaliana") {
    genmap[chr=="NC_003070.9", chr:=1]
    genmap[chr=="NC_003071.7", chr:=2]
    genmap[chr=="NC_003074.8", chr:=3]
    genmap[chr=="NC_003075.7", chr:=4]
    genmap[chr=="NC_003076.8", chr:=5]
    }
    if(species[[sppNum]]$spp%like%"Crotalus") {
    genmap[,chr:=gsub("chr", "", chr)]
    }
    nlp[,start:=pos]
    nlp[,stop:=pos+1]

    setkey(genmap, chr, start, stop)
    setkey(nlp, chr, start, stop)


    nlp.genmap <- foverlaps(nlp, genmap, nomatch=NA)

    dim(nlp)
    dim(nlp.genmap)

    nlp.genmap <- nlp.genmap[!duplicated(nlp.genmap, by="variant")]
    summary(nlp.genmap$V5)

    nlp <- nlp.genmap[,-c("start", "stop", "i.start", "i.stop")]
    setnames(nlp, c("V4", "V5"), c("genmap_id", "genmap_score"))

    message(paste("saving NLP + paramask/cnv + genmap dimensions:", dim(nlp), sep=" "))

    save(nlp, file=paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.paramask.genmap.Rdata", sep=""))
} else {
    message("loading genmap nlp")
    load(file=paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.paramask.genmap.Rdata", sep=""))
}
}

#############
### BUSCO ###
#############
if(!file.exists(paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.paramask.genmap.busco.Rdata", sep=""))) {
busco <- fread(species[[sppNum]]$BUSCO, fill=T, skip=3)
setnames(busco, c("V3", "V4", "V5", "V2"), c("chr", "start", "stop", "busco"))
nlp[,start:=pos]
nlp[,stop:=pos+1]

busco[V6=="-", stop_new:=start]
busco[V6=="-", start:=stop]
busco[V6=="-", stop:=stop_new]

if(species[[sppNum]]$spp%like%"Anopheles") {
    busco[chr=="NC_004818.2", chr:= "X" ]
    busco[chr=="NT_078265.2", chr:= "2L"]
    busco[chr=="NT_078266.2", chr:= "2R"]
    busco[chr=="NT_078267.5", chr:= "3L"]
    busco[chr=="NT_078268.4", chr:= "3R"]
}
if(species[[sppNum]]$spp%like%"Arabidopsis_thaliana") {
    busco[chr=="NC_003070.9", chr:=1]
    busco[chr=="NC_003071.7", chr:=2]
    busco[chr=="NC_003074.8", chr:=3]
    busco[chr=="NC_003075.7", chr:=4]
    busco[chr=="NC_003076.8", chr:=5]
}
if(species[[sppNum]]$spp%like%"Drosophila_simulans") {
    busco[chr=="NC_052520.2", chr:="sim_2L"]
    busco[chr=="NC_052521.2", chr:="sim_2R"]
    busco[chr=="NC_052522.2", chr:="sim_3L"]
    busco[chr=="NC_052523.2", chr:="sim_3R"]
    busco[chr=="NC_052524.2", chr:="sim_4" ]
    busco[chr=="NC_052525.2", chr:="sim_X" ]
}
if(species[[sppNum]]$spp%like%"Crotalus") {
    busco[,chr:=gsub("chr", "", chr)]
}

setkey(nlp, chr, start, stop)
setkey(busco, chr, start, stop)

nlp.busco <- foverlaps(nlp, busco[chr!=""][,c("chr", "start", "stop", "busco")], nomatch=NA)
nlp.busco <- nlp.busco[!duplicated(nlp.busco, by="variant")]

nlp <- nlp.busco[,-c("start", "stop", "i.start", "i.stop")]

message(paste("saving NLP + paramask/cnv + genmap + BUSCO dimensions:", dim(nlp), sep=" "))

save(nlp, file=paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.paramask.genmap.busco.Rdata", sep=""))

} else {
    message("loading busco nlp")
    load(file=paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.paramask.genmap.busco.Rdata", sep=""))
}

###
##########################
### make NLP aggregate ###
##########################
### version
agg_table_build_date <- "17_03_2026"  ### initial attemp; DMY
agg_table_nameS <- c("original", "posFIS", "noDev", "singleCopy", "multiCopy", "genmap", "busco", "no_adjacent_aa", "strict")

### patch
try(setnames(nlp, "col.x", "col", skip_absent=T))

### iterate
nlp.ag <- foreach(table_name=agg_table_nameS, .combine="rbind", .errorhandling="remove")%dopar%{

    if(table_name=="original") {
    filter <- rep(T, dim(nlp)[1])
    } else if(table_name=="posFIS") {
    filter <- nlp$FIS>0
    } else if(table_name=="noDev") {
    filter <- nlp$pval<=.05
    } else if(table_name=="singleCopy") {
    filter <- nlp$finalClass==0
    } else if(table_name=="multiCopy") {
    filter <- nlp$finalClass==1
    } else if(table_name=="genmap") {
    filter <- nlp$genmap_score==1
    } else if(table_name=="busco") {
    filter <- nlp$busco=="Complete"
    } else if(table_name=="no_adjacent_aa") {
    filter <- nlp$mutations_per_codon==1
    } else if(table_name=="strict") {
    filter <- nlp$mutations_per_codon==1 & nlp$busco=="Complete" & nlp$genmap_score==1 & nlp$finalClass==1
    }

    nlp.ag <- nlp[filter,list(nNS=sum(col=="missense_variant"), nS=sum(col=="synonymous_variant"), lNS=sum(nonsynonymous), lS=sum(synonymous)), list(nLocales_poly=nLocales_poly)]

    nlp.ag <- nlp.ag[order(nLocales_poly)]
    nlp.ag[,prop:=(nNS)/(nNS+nS)]
    nlp.ag[,spp:=spp]
    nlp.ag[,agg_build_date:=agg_table_build_date]
    nlp.ag[,agg_table_name:=table_name]
    nlp.ag[,NLP_table_build_date:=NLP_table_build_date]
    nlp.ag[,SNP_table_build_date:=SNP_table_build_date]
    nlp.ag[,prop:=nNS/(nNS+nS)]
    return(nlp.ag)
}
nlp.ag
nlp.ag[agg_table_name=="original"]

### output
write.csv(nlp.ag, file=paste(output_directory, spp, ".", agg_table_build_date, ".nlpTable_aggregate.csv", sep=""),
            quote=F, row.names=F)

seqClose(genofile)

