################################################################################################
# for the classification of same site, same nt, same AND different amino acids 
################################################################################################

library(data.table)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/all_variants_clean.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")

csv_class <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/classification/classificationScheme.csv")
class_dt <- read.csv(csv_class)

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
class_dt[class_dt$Classification == "A", "Count"] <- nrow(shared_classA)

### CLASS B: same position, same nt, diff amino acid ###########################################
shared_dt[
    ref_mel == ref_sim &
    alt_mel == alt_sim & 
    aa_ref_mel == aa_ref_sim & 
    aa_alt_mel != aa_alt_sim,
    classification := "B"
]
shared_classB <- shared_dt[classification == "B"]
message(nrow(shared_classB), " shared variants with same pos, nt, and diff aa")
class_dt[class_dt$Classification == "B", "Count"] <- nrow(shared_classB)

### CLASS C: same position, diff nt, same amino acid ############################################
shared_dt[
    ref_mel == ref_sim &
    alt_mel != alt_sim & 
    aa_ref_mel == aa_ref_sim & 
    aa_alt_mel == aa_alt_sim,
    classification := "C"
]
shared_classC <- shared_dt[classification == "C"]
message(nrow(shared_classC), " shared variants with same pos, diff nt, and same aa")
class_dt[class_dt$Classification == "C", "Count"] <- nrow(shared_classC)

### CLASS D: same position, diff nt, diff amino acid ###########################################
shared_dt[
    ref_mel == ref_sim &
    alt_mel != alt_sim & 
    aa_ref_mel == aa_ref_sim & 
    aa_alt_mel != aa_alt_sim,
    classification := "D"
]
shared_classD <- shared_dt[classification == "D"]
message(nrow(shared_classD), " shared variants with same pos, diff nt, and diff aa")
class_dt[class_dt$Classification == "D", "Count"] <- nrow(shared_classD)


############ UPDATE FILES ###############################
# outputs summary of classification stats 
table(shared_dt$classification, useNA = "ifany")

# RDS FILE
saveRDS(shared_dt, rds_file)

# CSV SUBSET
subset_table <- shared_dt[1:500, ]
table(subset_table$classification, useNA = "ifany")
fwrite(subset_table, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/all_variants_clean_500test.csv")
# CSV CLASSICATION
fwrite(class_dt, csv_class)

#########################################################


# ## CLASS X: testing new combos at same position ############################################
# shared_dt[
#     ref_mel == ref_sim &
#     alt_mel == alt_sim & 
#     aa_ref_mel == aa_ref_sim & 
#     aa_alt_mel != aa_alt_sim,
#     classification := "B"
# ]
# shared_dt[classification == "X", classification := NA]

# # shared_classX <- shared_dt[classification == "X"] # clear classification with class X
# message(nrow(shared_dt[classification == "X"]), " shared variants in class X")

# class_dt[class_dt$Classification=="H", "Count" ] <- nrow(shared_dt[classification == "H"])
