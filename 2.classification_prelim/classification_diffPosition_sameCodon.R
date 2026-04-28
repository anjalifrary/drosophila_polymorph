################################################################################################
# for the classification of different site but same codon, same / different nucleotides, same / different amino acids 
################################################################################################

library(data.table)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_variantsOfInterest_qualFiltered.rds")
test_rds_out <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_variantsOfInterest_qualFiltered_testing.rds")
shared_dt <- readRDS(test_rds_out)
# shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")

csv_class <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/classification/classificationScheme.csv")
class_dt <- fread(csv_class)
csv_test <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/classification/classificationScheme_testing.csv")


shared_dt[
    !is.na(ref_mel) & # start with mel rows, don't need to do sim bc all should reference each other but checking... 
    (
        adjacent_var_pos1 %in% shared_dt[!is.na(ref_sim), pos] |
        adjacent_var_pos2 %in% shared_dt[!is.na(ref_sim), pos] |
        adjacent_var_pos3 %in% shared_dt[!is.na(ref_sim), pos]
    ),
    classification_sameCodon := "I"
]
table(shared_dt$classification_sameCodon)

shared_dt[
    !is.na(ref_sim) & # checking sim, should be same
    (
        adjacent_var_pos1 %in% shared_dt[!is.na(ref_mel), pos] |
        adjacent_var_pos2 %in% shared_dt[!is.na(ref_mel), pos] |
        adjacent_var_pos3 %in% shared_dt[!is.na(ref_mel), pos]
    ),
    classification_sameCodon := "I"
]

table(shared_dt$classification_sameCodon) # should be same as before

# shared_dt[, classification_sameCodon := NA_character_]

shared_dt[, codon_id := paste(chr, codon_start_pos, sep = "_")]
unique_codons_I <- unique(
    shared_dt[classification_sameCodon == "I", .(chr, codon_start_pos)]
)
message(nrow(unique_codons_I))

class_dt <- rbind(
    class_dt,
    data.table(Classification = "I", Count = nrow(unique_codons_I)),
    fill = TRUE
)

###############################################################
### Analyzing the Data: Get number polymorphic sites within codons
###############################################################
# this gives a matrix where rows are n_mel meaning 0,1,2,or 3 polymorphic sites within codon
#  and columns are n_sim, similarly

# we are most interested in the class where n_mel and n_sim both = 1, so each species is polymorphic at 1 site within codon
codon_counts <- shared_dt[,
    .(n_mel = sum(!is.na(ref_mel)), n_sim = sum(!is.na(ref_sim))),
    by = .(chr, codon_start_pos)
]

joint_counts <- codon_counts[, .N, by = .(n_mel, n_sim)]

# idk how this works lol
joint_matrix <- dcast(
    joint_counts,
    n_mel ~ n_sim,
    value.var = "N",
    fill = 0
)

all_combos <- CJ(n_mel = 0:3, n_sim = 0:3)

joint_counts_full <- merge(
    all_combos,
    joint_counts,
    by = c("n_mel", "n_sim"),
    all.x = TRUE
)

joint_counts_full[is.na(N), N := 0]

joint_matrix <- dcast(
    joint_counts_full,
    n_mel ~ n_sim,
    value.var = "N"
)

setnames(
    joint_matrix,
    old = c("0", "1", "2", "3"),
    new = c("n_sim_0", "n_sim_1", "n_sim_2", "n_sim_3")
)
fwrite(joint_matrix,"/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/mel_sim_codonMatrix_NPolymorphicSitesWithinCodon.csv")


###############################################################
### Analyzing the Data: Get target codons where both species are polymorphic at ONE site within codon
###############################################################
 
# # get the codons
# target_codons <- codon_counts[
#     n_mel == 1 & n_sim == 1,
#     .(chr, codon_start_pos)
# ]

# # get the rows that correspond to variants within target codons
# target_dt <- shared_dt[
#     target_codons,
#     on = .(chr, codon_start_pos)
# ] # this should give 2 rows per codon 

target_dt <- shared_dt[
    ,
    if (
        .N == 2 &&
        sum(!is.na(ref_mel)) == 1 &&
        sum(!is.na(ref_sim)) == 1 &&
        uniqueN(pos) == 2   # ensures different positions (not same-site)
    ) .SD,
    by = .(chr, codon_start_pos)
]

target_dt[, .N, by = .(chr, codon_start_pos)][, table(N)]

mel_rows <- target_dt[!is.na(ref_mel)]
sim_rows <- target_dt[!is.na(ref_sim)]

paired_dt <- target_dt[,
    { 
        mel_row <- .SD[!is.na(ref_mel)]
        sim_row <- .SD[!is.na(ref_sim)]
        .(
            nt_pos_mel = mel_row$pos,
            nt_pos_sim = sim_row$pos,

            codon_ref_mel = mel_row$codon_ref_mel,
            codon_alt_mel = mel_row$codon_alt_mel,
            codon_ref_sim = sim_row$codon_ref_sim,
            codon_alt_sim = sim_row$codon_alt_sim,

            aa_ref_mel = mel_row$aa_ref_mel,
            aa_alt_mel = mel_row$aa_alt_mel,
            aa_ref_sim = sim_row$aa_ref_sim,
            aa_alt_sim = sim_row$aa_alt_sim
        )
    }, by = .(chr, codon_start_pos)
]

# make codon all uppercase
paired_dt[, `:=`(
    codon_ref_mel = toupper(codon_ref_mel),
    codon_alt_mel = toupper(codon_alt_mel),
    codon_ref_sim = toupper(codon_ref_sim),
    codon_alt_sim = toupper(codon_alt_sim)
)]



nrow(paired_dt) # number of 1-1 codons where genomic position is different (2 rows)

###### CLASSIFICATION! #############################################

### The codons ###
paired_dt[, same_ref_codon := codon_ref_mel == codon_ref_sim] # is reference codon the same

paired_dt[, same_alt_codon := codon_alt_mel == codon_alt_sim] # is alternate codon the same

paired_dt[, swapped_ref_alt_codon := (
    ( (codon_ref_mel == codon_alt_sim) |
    (codon_alt_mel == codon_ref_sim)) 
    & (codon_ref_mel != codon_alt_mel) 
    & (codon_ref_sim != codon_alt_sim)
)] # 7.5k codons out of 39k


### The amino acids ###
paired_dt[, same_ref_aa := aa_ref_mel == aa_ref_sim] 
paired_dt[, same_alt_aa := aa_alt_mel == aa_alt_sim]

paired_dt[, swapped_ref_alt_aa := ((
    (aa_ref_mel == aa_alt_sim) |
    (aa_alt_mel == aa_ref_sim))
    & (aa_ref_mel != aa_alt_mel) 
    & (aa_ref_sim != aa_alt_sim)
)]


message("-----STATS----- ")
table(paired_dt$same_ref_codon, useNA = "ifany")
table(paired_dt$same_alt_codon, useNA = "ifany")
table(paired_dt$swapped_ref_alt_codon, useNA = "ifany")

table(paired_dt$same_ref_aa, useNA = "ifany")
table(paired_dt$same_alt_aa, useNA = "ifany")
table(paired_dt$swapped_ref_alt_aa, useNA = "ifany")

names(paired_dt)

paired_dt[
    same_ref_codon == TRUE &
    same_alt_codon == FALSE &
    same_ref_aa == TRUE &
    same_alt_aa == TRUE,
    classification_sameCodon := "G"
]
shared_classX <- paired_dt[classification_sameCodon == "G"]
message(nrow(shared_classX), " shared variants with diff pos, same codon, and aa")
class_dt$Count[class_dt$Classification == "G"] <- nrow(shared_classX)
shared_classX <- NULL
table(paired_dt$classification_sameCodon, useNA = "ifany")

###

paired_dt[
    same_ref_codon == FALSE &
    same_alt_codon == TRUE &
    same_ref_aa == TRUE &
    same_alt_aa == TRUE,
    classification_sameCodon := "H"
]
shared_classX <- paired_dt[classification_sameCodon == "H"]
message(nrow(shared_classX), " shared variants with diff pos, same codon, and aa")
class_dt$Count[class_dt$Classification == "H"] <- nrow(shared_classX)
shared_classX <- NULL
table(paired_dt$classification_sameCodon, useNA = "ifany")

###

paired_dt[
    same_ref_codon == TRUE &
    same_alt_codon == FALSE &
    same_ref_aa == TRUE &
    same_alt_aa == FALSE,
    classification_sameCodon := "I"
]
shared_classX <- paired_dt[classification_sameCodon == "I"]
message(nrow(shared_classX), " shared variants with diff pos, same codon, and aa")
class_dt$Count[class_dt$Classification == "I"] <- nrow(shared_classX)
shared_classX <- NULL
table(paired_dt$classification_sameCodon, useNA = "ifany")

### 

paired_dt[
    same_ref_codon == FALSE &
    same_alt_codon == TRUE &
    same_ref_aa == FALSE &
    same_alt_aa == TRUE,
    classification_sameCodon := "J"
]
shared_classX <- paired_dt[classification_sameCodon == "J"]
message(nrow(shared_classX), " shared variants with diff pos, same codon, and aa")
class_dt$Count[class_dt$Classification == "J"] <- nrow(shared_classX)
shared_classX <- NULL
table(paired_dt$classification_sameCodon, useNA = "ifany")

###

paired_dt[
    same_ref_codon == FALSE &
    same_alt_codon == FALSE &
    same_ref_aa == TRUE &
    same_alt_aa == TRUE,
    classification_sameCodon := "K"
]
shared_classX <- paired_dt[classification_sameCodon == "K"]
message(nrow(shared_classX), " shared variants with diff pos, same codon, and aa")
class_dt$Count[class_dt$Classification == "K"] <- nrow(shared_classX)
shared_classX <- NULL
table(paired_dt$classification_sameCodon, useNA = "ifany")

###

paired_dt[
    same_ref_codon == FALSE &
    same_alt_codon == FALSE &
    same_ref_aa == TRUE &
    same_alt_aa == FALSE,
    classification_sameCodon := "L"
]
shared_classX <- paired_dt[classification_sameCodon == "L"]
message(nrow(shared_classX), " shared variants with diff pos, same codon, and aa")
class_dt$Count[class_dt$Classification == "L"] <- nrow(shared_classX)
shared_classX <- NULL
table(paired_dt$classification_sameCodon, useNA = "ifany")

###

paired_dt[
    same_ref_codon == FALSE &
    same_alt_codon == FALSE &
    same_ref_aa == FALSE &
    same_alt_aa == TRUE,
    classification_sameCodon := "M"
]

shared_classX <- paired_dt[classification_sameCodon == "M"]
message(nrow(shared_classX), " shared variants with diff pos, same codon, and aa")
class_dt$Count[class_dt$Classification == "M"] <- nrow(shared_classX)
shared_classX <- NULL
table(paired_dt$classification_sameCodon, useNA = "ifany")
###

paired_dt[
    same_ref_codon == FALSE &
    same_alt_codon == FALSE &
    same_ref_aa == FALSE &
    same_alt_aa == FALSE,
    classification_sameCodon := "N"
]
shared_classX <- paired_dt[classification_sameCodon == "N"]
message(nrow(shared_classX), " shared variants with diff pos, same codon, and aa")
class_dt$Count[class_dt$Classification == "N"] <- nrow(shared_classX)
shared_classX <- NULL
table(paired_dt$classification_sameCodon, useNA = "ifany")

### for the swapped codons / amino acids:

# EXCLUSIVE OR 
# changing same_ref_codon and same_alt_codon ONLY for a-d suffices 
paired_dt[
    (xor(
        (codon_ref_mel == codon_alt_sim),
        (codon_alt_mel == codon_ref_sim)
    )) &
    same_ref_aa == FALSE & 
    same_alt_aa == FALSE,
    classification_sameCodon := ifelse(
        is.na(classification_sameCodon),
        "d",
        paste0(classification_sameCodon, ";d")
    )
]

table(paired_dt$classification_sameCodon, useNA = "ifany")

### 

# for the full codon SWAP - NONE
paired_dt[
    ((codon_ref_mel == codon_alt_sim) & (codon_alt_mel == codon_ref_sim)) &
    same_ref_aa == FALSE & 
    same_alt_aa == TRUE,
    classification_sameCodon := ifelse(
        is.na(classification_sameCodon),
        "e",
        paste0(classification_sameCodon, ";e")
    )
]

table(paired_dt$classification_sameCodon, useNA = "ifany")


# EXCLUSIVE OR 
# changing same_ref_aa and same_alt_aa ONLY for e-h suffices 
paired_dt[
    (xor(
        (aa_ref_mel == aa_alt_sim),
        (aa_alt_mel == aa_ref_sim)
    )) &
    same_ref_codon == FALSE & 
    same_alt_codon == TRUE,
    classification_sameCodon := ifelse(
        is.na(classification_sameCodon),
        "g",
        paste0(classification_sameCodon, ";g")
    )]
table(paired_dt$classification_sameCodon, useNA = "ifany")

# bi-directional swap

paired_dt[
    ((aa_ref_mel == aa_alt_sim) & (aa_alt_mel == aa_ref_sim)) &
    same_ref_codon == TRUE & 
    same_alt_codon == TRUE,
    classification_sameCodon := ifelse(
        is.na(classification_sameCodon),
        "k",
        paste0(classification_sameCodon, ";k")
    )]
table(paired_dt$classification_sameCodon, useNA = "ifany")

# testing both codon and aa are swapped
paired_dt[
    ((aa_ref_mel == aa_alt_sim) & (aa_alt_mel == aa_ref_sim)) &
    swapped_ref_alt_codon == TRUE,
    classification_sameCodon := ifelse(
        is.na(classification_sameCodon),
        "k",
        paste0(classification_sameCodon, ";k")
    )]
table(paired_dt$classification_sameCodon, useNA = "ifany")

#     G   G;h   H;i     I   I;f     J   J;g      K K;a;j   K;j     L L;b;e   L;e
#    36  1736    50  5570 21832    43   152     1   593    46   313  5067   232
#     M M;c;e   M;e     N N;d;e N;d;j   N;e   N;j
#    14   847    38  1112  1036     8    75     1

#   G     G;h     H;i       I     I;f       J     J;g       K K;a;j;k     K;j
#  36    1736      50    5570   21832      43     152       1     593      46
#   L   L;b;e     L;e       M   M;c;e     M;e       N   N;d;e N;d;j;k     N;e
# 313    5067     232      14     847      38    1112    1036       8      75
# N;j
#   1

# G     = 3 codons, 2 AA
# G;h   = 3 codons, 1 AA
# H;i   = 3 codons, 1 AA
# I     = 3 codons, 3 AA
# I;f   = 3 codons, 2 AA
# J     = 3 codons, 3 AA
# J;g   = 3 codons, 2 AA
# K     = 4 codons, 2 AA
# K;a;j = 3 codons, 1 AA
# K;j   = 4 codons, 1 AA
# L     = 4 codons, 3 AA
# L;b;e = 3 codons, 2 AA 
# L;e   = 4 codons, 2 AA
# M     = 4 codons, 3 AA
# M;c;e = 3 codons, 2 AA
# M;e   = 4 codons, 2 AA
# N     = 4 codons, 4 AA
# N;d;e = 3 codons, 3 AA
# N;d;j = 3 codons, 2 AA
# N;e   = 4 codons, 3 AA
# N;j   = 4 codons, 2 AA

############ UPDATE FILES ###############################

######## test files!!! ##################
# outputs summary of classification stats 
table(shared_dt$classification_sameCodon, useNA = "ifany")
table(shared_dt$classification_sameSite, useNA = "ifany")

shared_dt[, .N, by = .(classification_sameSite, classification_sameCodon)]

# RDS FILE
saveRDS(shared_dt, test_rds_out)

# CSV SUBSET
subset_table <- shared_dt[1:500, ]
table(subset_table$classification, useNA = "ifany")
fwrite(subset_table, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/subset_variantsOfInterest_qualFiltered_500test.csv")
# CSV CLASSICATION
fwrite(class_dt, csv_test)
