################################################################################################
# for the classification of same site, same / different nucleotides, same / different amino acids 
################################################################################################

library(data.table)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_variantsOfInterest_qualFiltered.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")

csv_class <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/classification/classificationScheme.csv")
class_dt <- fread(csv_class)

### CLASS A: same position, same nt, same amino acid ############################################
shared_dt[
    ref_mel == ref_sim &
    alt_mel == alt_sim & 
    aa_ref_mel == aa_ref_sim & 
    aa_alt_mel == aa_alt_sim,
    classification := "A"
]
shared_classA <- shared_dt[classification == "A"]
message(nrow(shared_classA), " shared variants with same pos, nt, and aa")
class_dt$Count[class_dt$Classification == "A"] <- nrow(shared_classA)

### CLASS B: same position, same nt, diff amino acid ###########################################
shared_dt[
    ref_mel == ref_sim &
    alt_mel == alt_sim & 
    aa_ref_mel != aa_ref_sim & 
    aa_alt_mel == aa_alt_sim,
    classification := "B2"
]
shared_classB <- shared_dt[classification == "B2"]
message(nrow(shared_classB), " shared variants with same pos, nt, and diff aa")
class_dt$Count[class_dt$Classification == "B2"] <- nrow(shared_classB)

shared_classB <- NULL
### CLASS C: same position, diff nt, same amino acid ############################################
shared_dt[
    ref_mel == ref_sim &
    alt_mel != alt_sim & 
    aa_ref_mel == aa_ref_sim & 
    aa_alt_mel == aa_alt_sim,
    classification := "C1"
]
# shared_classC <- NULL
shared_classC <- shared_dt[classification == "C1"]
message(nrow(shared_classC), " shared variants with same pos, diff nt, and same aa")
class_dt$Count[class_dt$Classification == "C1"] <- nrow(shared_classC)

### CLASS D: same position, diff nt, diff amino acid ###########################################
shared_dt[
    ref_mel == ref_sim &
    alt_mel != alt_sim & 
    aa_ref_mel == aa_ref_sim & 
    aa_alt_mel != aa_alt_sim,
    classification := "D1"
]
shared_classD <- shared_dt[classification == "D1"]
message(nrow(shared_classD), " shared variants with same pos, diff nt, and diff aa")
class_dt$Count[class_dt$Classification == "D1"] <- nrow(shared_classD)
shared_classD <- NULL


## CLASS X: testing new combos at same position ############################################
shared_dt[
    ref_mel == ref_sim &
    alt_mel == alt_sim & 
    aa_ref_mel != aa_ref_sim & 
    aa_alt_mel != aa_alt_sim,
    classification := "E"
]
shared_classX <- shared_dt[classification == "E"]
message(nrow(shared_classX), " shared variants in class X")
class_dt$Count[class_dt$Classification=="E"] <- nrow(shared_classX)
shared_classX <- NULL

############ UPDATE FILES ###############################
# outputs summary of classification stats 
table(shared_dt$classification_sameSite, useNA = "ifany")

# RDS FILE
saveRDS(shared_dt, rds_file)

# CSV SUBSET
subset_table <- shared_dt[1:500, ]
table(subset_table$classification, useNA = "ifany")
fwrite(subset_table, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variantsOfInterest_qualFiltered_500test.csv")
# CSV CLASSICATION
fwrite(class_dt, csv_class)

#########################################################

