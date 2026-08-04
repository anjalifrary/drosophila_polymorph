# library(SeqArray)
library(data.table)
# library(doMC)
# registerDoMC(16)

results <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/adaptedMK/asymptotic_MK_longResults_polyAF_speciesSpecificBG.rds")

nrow(results[!is.na(alpha)])

