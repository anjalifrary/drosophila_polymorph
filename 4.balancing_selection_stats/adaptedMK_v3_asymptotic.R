library(SeqArray)
library(data.table)
library(doMC)
registerDoMC(16)


# bg_rds <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/all_quality_variants_MAF5_merge_unfilt.rds")
nrow(bg_rds)

# classed_rds <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/voi_fromBG_qualVar_ofInterest_MAF5_classed_06-29-2026.rds") # only classes ABFGOPXY
nrow(classed_rds)
names(classed_rds)

# candidate files of chr + pos (txt) at /scratch/ejy4bu/drosophila/GO/gowinda/candidateFiles/MAFxxfilter_polyAF/candidate_chrpos_AB_xx_polyAF.txt
# bg files of chr + pos (.txt) at /scratch/ejy4bu/drosophila/GO/gowinda/backgroundFiles/MAFxxfilter_polyAF/bg_speciesSpecific_xx_polyAF.txt

tsp <- c("A", "B")
conv <- c("F", "G", "O", "P", "X", "Y")


# use bg and candidate files like the ones made for gowinda (but don't save as text files)
# add gene info for per-gene stats 

# run at varying MAF 
# how to save it? 