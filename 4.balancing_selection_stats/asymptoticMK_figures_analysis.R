# library(SeqArray)
library(data.table)
# library(doMC)
# registerDoMC(16)

results <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/adaptedMK/asymptotic_MK_longResults_polyAF_speciesSpecificBG.rds")

nrow(results[!is.na(alpha)])

summary(results[, alpha])
summary(results[, logOR])

### FIGURES ###

### 1. histogram, breaks=100 
# hist symmetry should shift to center around 0 as MAF increases... ? 
# add -inf / inf bar on the ends ? 
maf = 0.005
group = "TSP"
MAF_def = "polyAF"


### 2. box-whisker plots, alpha vs MAF with annotated label# genes for which alpha !is.na

### 3. line graph alpha vs MAF where each line is a gene (like asymptotic application)

### 4. line graph log(OR) vs MAF where each line is a gene (like asymptotic application)

