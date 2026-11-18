### FROM ALAN - implemented these updates into poolseq_nlp_polyAF_11.2026.R


# ijob -A berglandlab_standard -c1 -p standard --mem=20G
# module load gcc/14.2.0 openmpi/5.0.7 R/4.6.0; R


### which species to work on
  args <- commandArgs(trailingOnly = TRUE)
  sppNum <- as.numeric(args[1])

  # sppNum <- 7

### define directory
  output_directory <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/nlp/"

### output file stems
  NLP_table_build_date <- "16_11_2026"  ### D_M_Y
  agg_table_build_date <- "16_11_2026"  ### D_M_Y

### libraries
  library(SeqArray)
  library(data.table)
  #library(ggplot2)
  library(foreach)
  library(doMC)
  registerDoMC(10)
  library(readxl)
  library(SNPRelate)
  #library(poolfstat)
  # 103

### generate meta-meta object
#   source("~/misc/1000G/done/1.collectiveAnalysis/0.createMetadata.R")
#   names(species)[sppNum]
#   spp <- species[[sppNum]]$spp

### open GDS file
  try(genofile <- seqOpen(species[[sppNum]]$gdsFile, allow.duplicate=TRUE))

### generate metadata file
  source(species[[sppNum]]$sample_metadata)
  str(samps)

### what species?
  message(paste("working on species: ", spp))

######################
### make SNP table ###
######################
  ### version
    SNP_table_build_date <- "16_11_2026"  ### initial attemp; DMY

  ### run or load
    if(!file.exists(paste(output_directory, spp, ".", SNP_table_build_date, ".snpTable.Rdata", sep=""))) {
      source("~/misc/1000G/done/1.collectiveAnalysis/1.makeSNP_table.R")
      snp.dt[,spp:=spp]
      save(snp.dt, file=paste(output_directory, spp, ".", SNP_table_build_date, ".snpTable.Rdata", sep=""))
    } else {
      message("Current SNP table already exists; loading")
      load(paste(output_directory, spp, ".", SNP_table_build_date, ".snpTable.Rdata", sep=""))
    }

  ### inspect
    str(snp.dt)
    setkey(snp.dt, col)

######################
### make NLP table ###
######################
  ### version
  
  ### run or load
    source("~/misc/1000G/done/1.collectiveAnalysis/2.makeNLP_table.R")

    if(!file.exists(paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.Rdata", sep=""))) {

      nlp <- foreach(bin.i=unique(snp.dt$bin), .errorhandling="remove")%dopar%{
        # bin.i <- 7
        nlp_fun(bin.i=bin.i, lib_type=species[[sppNum]]$libraryType, colTypes=species[[sppNum]]$snpType)
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
        cnv <- as.data.table(read_excel(paste(github_root, github_stem, "arabidopsis/TPC2019-LSB-00640R2_Supplemental_Data_Set_1_8.xlsx", sep="/"), skip=4, sheet="Supplemental Data Set 1"))
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

      if(species[[sppNum]]$spp=="Mimulus_gutattus") {
        genmap[chr%like%"Chr", chr:=as.numeric(tstrsplit(chr, "_")[[2]])]
        table(genmap$chr)
      }
      if(species[[sppNum]]$spp=="Drosophila_suzukii") {
        genmap[,chr:=gsub("chr", "", chr)]
      }
      if(species[[sppNum]]$spp%like%"Zaprionus") {
        genmap[chr=="JAUIZU010000001.1", chr:="Scaffold_1"]
        genmap[chr=="JAUIZU010000002.1", chr:="Scaffold_2"]
        genmap[chr=="JAUIZU010000003.1", chr:="Scaffold_3"]
        genmap[chr=="JAUIZU010000004.1", chr:="Scaffold_4"]
        genmap[chr=="JAUIZU010000005.1", chr:="Scaffold_5"]
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
    if(species[[sppNum]]$spp=="Mimulus_gutattus") {
        busco[chr%like%"Chr", chr:=as.numeric(tstrsplit(chr, "_")[[2]])]
    }
    if(species[[sppNum]]$spp=="Drosophila_suzukii") {
        busco[,chr:=gsub("chr", "", chr)]
    }
    if(species[[sppNum]]$spp%like%"Zaprionus") {
        busco[chr=="JAUIZU010000001.1", chr:="Scaffold_1"]
        busco[chr=="JAUIZU010000002.1", chr:="Scaffold_2"]
        busco[chr=="JAUIZU010000003.1", chr:="Scaffold_3"]
        busco[chr=="JAUIZU010000004.1", chr:="Scaffold_4"]
        busco[chr=="JAUIZU010000005.1", chr:="Scaffold_5"]
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

#### add in repeat masker & wmDust annotations
  if(!file.exists(paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.paramask.genmap.busco.repmask.wd.Rdata", sep=""))) {
    repmask <- fread(species[[sppNum]]$rm, fill=T, skip=2)
    repmask <- repmask[,c("V5", "V6", "V7", "V11")]
    setnames(repmask, names(repmask), c("chr", "start", "stop", "repmask_call"))
    repmask[stop<start, newStop:=start]
    repmask[stop<start, newStart:=stop]
    repmask[!is.na(newStop), stop:=newStop]
    repmask[!is.na(newStart), start:=newStart]
    repmask <- repmask[,-c("newStart", "newStop")]
    repmask[stop<start]

    wmdust <- fread(species[[sppNum]]$wm, fill=T, skip=0)
    wmdust[,V11:="wmdust"]
    setnames(wmdust, names(wmdust), c("chr", "start", "stop", "wmdust_call"))
    wmdust[stop<start]
    wmdust[stop<start, newStop:=start]
    wmdust[stop<start, newStart:=stop]
    wmdust[!is.na(newStop), stop:=newStop]
    wmdust[!is.na(newStart), start:=newStart]
    wmdust <- wmdust[,-c("newStart", "newStop")]
    wmdust[stop<start]


     if(species[[sppNum]]$spp%like%"Anopheles") {
      repmask[chr=="NC_004818.2", chr:= "X" ]
      repmask[chr=="NT_078265.2", chr:= "2L"]
      repmask[chr=="NT_078266.2", chr:= "2R"]
      repmask[chr=="NT_078267.5", chr:= "3L"]
      repmask[chr=="NT_078268.4", chr:= "3R"]
    }
    if(species[[sppNum]]$spp%like%"Arabidopsis_thaliana") {
      repmask[chr=="NC_003070.9", chr:=1]
      repmask[chr=="NC_003071.7", chr:=2]
      repmask[chr=="NC_003074.8", chr:=3]
      repmask[chr=="NC_003075.7", chr:=4]
      repmask[chr=="NC_003076.8", chr:=5]
    }
    if(species[[sppNum]]$spp%like%"Drosophila_simulans") {
      repmask[chr=="NC_052520.2", chr:="sim_2L"]
      repmask[chr=="NC_052521.2", chr:="sim_2R"]
      repmask[chr=="NC_052522.2", chr:="sim_3L"]
      repmask[chr=="NC_052523.2", chr:="sim_3R"]
      repmask[chr=="NC_052524.2", chr:="sim_4" ]
      repmask[chr=="NC_052525.2", chr:="sim_X" ]
    }
    if(species[[sppNum]]$spp%like%"Crotalus") {
      repmask[,chr:=gsub("chr", "", chr)]
    }
    if(species[[sppNum]]$spp=="Mimulus_gutattus") {
        repmask[chr%like%"Chr", chr:=as.numeric(tstrsplit(chr, "_")[[2]])]
    }
    if(species[[sppNum]]$spp=="Drosophila_suzukii") {
        repmask[,chr:=gsub("chr", "", chr)]
    }
    if(species[[sppNum]]$spp%like%"Zaprionus") {
        repmask[chr=="JAUIZU010000001.1", chr:="Scaffold_1"]
        repmask[chr=="JAUIZU010000002.1", chr:="Scaffold_2"]
        repmask[chr=="JAUIZU010000003.1", chr:="Scaffold_3"]
        repmask[chr=="JAUIZU010000004.1", chr:="Scaffold_4"]
        repmask[chr=="JAUIZU010000005.1", chr:="Scaffold_5"]
    }


    if(species[[sppNum]]$spp%like%"Anopheles") {
      wmdust[chr=="NC_004818.2", chr:= "X" ]
      wmdust[chr=="NT_078265.2", chr:= "2L"]
      wmdust[chr=="NT_078266.2", chr:= "2R"]
      wmdust[chr=="NT_078267.5", chr:= "3L"]
      wmdust[chr=="NT_078268.4", chr:= "3R"]
    }
    if(species[[sppNum]]$spp%like%"Arabidopsis_thaliana") {
      wmdust[chr=="NC_003070.9", chr:=1]
      wmdust[chr=="NC_003071.7", chr:=2]
      wmdust[chr=="NC_003074.8", chr:=3]
      wmdust[chr=="NC_003075.7", chr:=4]
      wmdust[chr=="NC_003076.8", chr:=5]
    }
    if(species[[sppNum]]$spp%like%"Drosophila_simulans") {
      wmdust[chr=="NC_052520.2", chr:="sim_2L"]
      wmdust[chr=="NC_052521.2", chr:="sim_2R"]
      wmdust[chr=="NC_052522.2", chr:="sim_3L"]
      wmdust[chr=="NC_052523.2", chr:="sim_3R"]
      wmdust[chr=="NC_052524.2", chr:="sim_4" ]
      wmdust[chr=="NC_052525.2", chr:="sim_X" ]
    }
    if(species[[sppNum]]$spp%like%"Crotalus") {
      wmdust[,chr:=gsub("chr", "", chr)]
    }
    if(species[[sppNum]]$spp=="Mimulus_gutattus") {
        wmdust[chr%like%"Chr", chr:=as.numeric(tstrsplit(chr, "_")[[2]])]
    }
    if(species[[sppNum]]$spp=="Drosophila_suzukii") {
        wmdust[,chr:=gsub("chr", "", chr)]
    }
    if(species[[sppNum]]$spp%like%"Zaprionus") {
        wmdust[chr=="JAUIZU010000001.1", chr:="Scaffold_1"]
        wmdust[chr=="JAUIZU010000002.1", chr:="Scaffold_2"]
        wmdust[chr=="JAUIZU010000003.1", chr:="Scaffold_3"]
        wmdust[chr=="JAUIZU010000004.1", chr:="Scaffold_4"]
        wmdust[chr=="JAUIZU010000005.1", chr:="Scaffold_5"]
    }
    nlp[,start:=pos]
    nlp[,stop:=pos+1]

    setkey(nlp, chr, start, stop)
    setkey(repmask, chr, start, stop)
    setkey(wmdust, chr, start, stop)

    nlp.repmask <- foverlaps(nlp, repmask[chr!=""][,c("chr", "start", "stop", "repmask_call")], nomatch=NA)
    nlp.repmask <- nlp.repmask[!duplicated(nlp.repmask, by="variant")]

    nlp <- nlp.repmask[,-c("start", "stop", "i.start", "i.stop")]
    nlp[,start:=pos]
    nlp[,stop:=pos+1]
    setkey(nlp, chr, start, stop)

    nlp.wmdust <- foverlaps(nlp, wmdust[chr!=""][,c("chr", "start", "stop", "wmdust_call")], nomatch=NA)
    nlp.wmdust <- nlp.wmdust[!duplicated(nlp.wmdust, by="variant")]

    nlp <- nlp.wmdust[,-c("start", "stop", "i.start", "i.stop")]
    
    message(paste("saving NLP + paramask/cnv + genmap + BUSCO + repmask + wmdust dimensions:", dim(nlp), sep=" "))

    save(nlp, file=paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.paramask.genmap.busco.repeatmask.wmdust.Rdata", sep=""))
  } else {
      message("loading repeatmask/wmdust nlp")
      load(file=paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.paramask.genmap.busco.repeatmask.wmdust.Rdata", sep=""))
  }


### make down-projection
    message("down-projection")
    nlp.miss <- nlp[,list(.N), nLocales_total]
    nlp.miss <- nlp.miss[order(-nLocales_total)]
    nlp.miss[,csN:=cumsum(N)]
    nlp.miss[,fcsN:=csN/max(nlp.miss$csN)]
    
    nlp.miss
    projTarget <- nlp.miss[which.min(abs(fcsN-0.95))]$nLocales_total

    nlp[nLocales_total>=projTarget, nLocales_poly_proj_round:=round((nLocales_poly/nLocales_total) * projTarget)]
    nlp[nLocales_total>=projTarget, nLocales_poly_proj_binomial:=rbinom(dim(nlp[nLocales_total>=projTarget])[1], projTarget, (nLocales_poly/nLocales_total))]

    table(is.na(nlp$nLocales_poly_proj_round), nlp$nLocales_total)
    save(nlp, file=paste(output_directory, spp, ".", NLP_table_build_date, ".nlpTable.paramask.genmap.busco.repeatmask.wmdust.projection.Rdata", sep=""))

### close up and finish
  message("done")
  seqClose(genofile)


###
###########################
#### make NLP aggregate ###
###########################
#  ### version
#    # agg_table_names <- c("original", "posFIS", "noDev", "singleCopy", "multiCopy", "genmap", "busco", "no_adjacent_aa", "strict")
#    agg_tables <- list()
#      agg_tables$unfiltered <- list()
#      agg_tables$unfiltered$filter <- rep(T, dim(nlp)[1])
#      agg_tables$unfiltered$name <- "unfiltered"
#
#    agg_tables$basic_filter <- list()
#    agg_tables$basic_filter$name <- "basic_filter"
#
#    agg_tables$basic_filter$filter <- is.na(nlp$wmdust_call) &
#                                      is.na(nlp$repmask_call) &
#                                      nlp$genmap_score==1 &
#                                      nlp$mutations_per_codon==1 & 
#                                      nlp$use==T
#
#
#    agg_tables$strict_filter <- list()
#    agg_tables$strict_filter$name <- "strict_filter"
#    agg_tables$strict_filter$filter <-  is.na(nlp$wmdust_call) &
#                                          is.na(nlp$repmask_call) &
#                                          nlp$genmap_score==1 &
#                                          nlp$mutations_per_codon==1 &
#                                          nlp$use==T &
#                                          nlp$busco=="Complete"
#
#
#  lapply(agg_tables, function(x) table(x[2]$filter))
#
#
#  ### patch
#    try(setnames(nlp, "col.x", "col", skip_absent=T))
#
#  ### iterate
#    nlp.ag <- foreach(filter_table=agg_tables, .combine="rbind", .errorhandling="remove")%dopar%{
#      # filter_table <- agg_tables[[1]]
#
#
#      nlp.ag <- nlp[filter_table$filter,list(nNS=sum(col=="missense_variant"), nS=sum(col=="synonymous_variant"), lNS=sum(nonsynonymous), lS=sum(synonymous)), list(nLocales_poly=nLocales_poly)]
#
#      nlp.ag <- nlp.ag[order(nLocales_poly)]
#      nlp.ag[,prop:=(nNS)/(nNS+nS)]
#      nlp.ag[,spp:=spp]
#      nlp.ag[,agg_build_date:=agg_table_build_date]
#      nlp.ag[,agg_table_name:= filter_table$name ]
#      nlp.ag[,NLP_table_build_date:=NLP_table_build_date]
#      nlp.ag[,SNP_table_build_date:=SNP_table_build_date]
#      nlp.ag[,prop:=nNS/(nNS+nS)]
#      return(nlp.ag)
#    }
#    nlp.ag
#    table(nlp.ag$agg_table_name)
#
#  ### output
#    write.csv(nlp.ag, file=paste(output_directory, spp, ".", agg_table_build_date, ".nlpTable_aggregate.csv", sep=""),
#              quote=F, row.names=F)

