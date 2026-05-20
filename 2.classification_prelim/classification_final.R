library(data.table)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_final.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total rows in shared table")

csv_class <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/classification/classification_table_final.csv")


# function to get an 'unordered' set for codon and amino acid comparison (ordering by alphabetization)
get_pair <- function(ref,alt) {
    if(is.na(ref) | is.na(alt)) return(NA_character_)
    sort(c(ref,alt)) # returns a vector
}

# classify <- function(mel, sim) {
#     if (any(is.na(mel)) | any(is.na(sim)))  return(NA_character_)
#     if (setequal(mel,sim))                  return("both_shared")
#     if (length(intersect(mel, sim)) == 1)   return("one_shared")
#     return("zero_shared")
# }

count_shared <- function(mel, sim){
    if (any(is.na(mel)) | any(is.na(sim))){
        return(c(
            total = NA_integer_,
            shared = NA_integer_,
            Nshared = NA_integer_
        ))
    }
    mel <- unique(mel)
    sim <- unique(sim)
    total <- length(union(mel,sim))
    shared <- length(intersect(mel,sim))
    c(
        total = total,
        shared = shared,
        Nshared = total-shared
    )
}

# get codon pairs for mel and sim: 
shared_dt[, codon_pair_mel := mapply(
    get_pair,
    codon_ref_mel,
    codon_alt_mel,
    SIMPLIFY = FALSE
)]
shared_dt[, aa_pair_mel := mapply(
    get_pair,
    aa_ref_mel,
    aa_alt_mel,
    SIMPLIFY = FALSE
)]

shared_dt[, codon_pair_sim := mapply(
    get_pair,
    codon_ref_sim,
    codon_alt_sim,
    SIMPLIFY = FALSE
)]
shared_dt[, aa_pair_sim := mapply(
    get_pair,
    aa_ref_sim,
    aa_alt_sim,
    SIMPLIFY = FALSE
)]

# codons:
shared_dt[, c("total_codons", "shared_codons","Nshared_codons") :=
    transpose(
        mapply(
            count_shared,
            codon_pair_mel,
            codon_pair_sim,
            SIMPLIFY = FALSE
        )
    )
]
# amino acids:
shared_dt[, c("total_aa", "shared_aa", "Nshared_aa") :=
    transpose(
        mapply(
            count_shared,
            aa_pair_mel,
            aa_pair_sim,
            SIMPLIFY = FALSE
        )
    )
]

shared_dt[, mel_aa := lengths(lapply(aa_pair_mel, unique))]
shared_dt[, sim_aa := lengths(lapply(aa_pair_sim, unique))]


#### same site
same_site <- shared_dt[
    !is.na(ref_mel) & !is.na(ref_sim)
]   
same_site[, same_pos := 1L]


same_summary <- same_site[
    ,
    .(
        chr,
        codon_start_pos,
        same_pos,
        total_codons,
        shared_codons,
        Nshared_codons,
        total_aa,
        shared_aa,
        Nshared_aa,
        mel_aa,
        sim_aa
    )
]



################### for different site variants:
diff_site <- shared_dt[
    (!is.na(ref_mel) & is.na(ref_sim)) |
    (!is.na(ref_sim) & is.na(ref_mel))
]
message(nrow(diff_site), "rows for diff-site variants")
# diff_site now contains collapsed table where each row is a codon; all rows were duplicated... 
# will then assign classification column and merge back to shared_dt on chr and codon_start_pos

# collapse to one row per codon
diff_summary <- diff_site[,
    {
        mel_row <- which(!is.na(ref_mel))[1]
        sim_row <- which(!is.na(ref_sim))[1]

        # codon counts
        codon_counts <- count_shared(
            codon_pair_mel[[mel_row]],
            codon_pair_sim[[sim_row]]
        )

        # amino acid counts
        aa_counts <- count_shared(
            aa_pair_mel[[mel_row]],
            aa_pair_sim[[sim_row]]
        )
        .(
            same_pos = 0L,
            total_codons   = codon_counts["total"],
            shared_codons  = codon_counts["shared"],
            Nshared_codons = codon_counts["Nshared"],

            total_aa   = aa_counts["total"],
            shared_aa  = aa_counts["shared"],
            Nshared_aa = aa_counts["Nshared"],

            mel_aa = length(unique(aa_pair_mel[[mel_row]])),
            sim_aa = length(unique(aa_pair_sim[[sim_row]]))
        )
    },
    by = .(chr, codon_start_pos)
]

# merge same and diff site variants
all_classes <- rbindlist(
    list(same_summary, diff_summary),
    use.names = TRUE
)



###### auto-build table

class_table <- all_classes[
    ,
    .(Count = .N),
    by = .(
        same_pos,
        total_codons,
        shared_codons,
        Nshared_codons,
        total_aa,
        shared_aa,
        Nshared_aa,
        mel_aa,
        sim_aa
    )
]

setorder(
    class_table,
    -same_pos,          # descending order
    total_codons,
    shared_codons,
    -total_aa,
    shared_aa
)

class_table[, Classification := LETTERS[seq_len(.N)]]

setcolorder(
    class_table,
    c(
        "Classification",
        "Count",
        "same_pos",
        "total_codons",
        "shared_codons",
        "Nshared_codons",
        "total_aa",
        "shared_aa",
        "Nshared_aa",
        "mel_aa",
        "sim_aa"
    )
)

# commented out to avoid overwriting existing annotated table
# # save class_table
# fwrite(class_table, csv_class)
# message("classification table written to: ", csv_class)


# update shared table classification column
class_cols <- c(
    "same_pos",
    "total_codons",
    "shared_codons",
    "Nshared_codons",
    "total_aa",
    "shared_aa",
    "Nshared_aa",
    "mel_aa",
    "sim_aa"
)

all_classes <- merge(
    all_classes,
    class_table[, c("Classification", class_cols), with = FALSE],
    by = class_cols,
    all.x = TRUE
)

# clear class column
shared_dt[, classification := NA_character_]

# merge classification on chr and codon start pos
shared_dt[
    all_classes[, .(
        chr,
        codon_start_pos,
        Classification
    )],
    on = .(chr, codon_start_pos),
    classification := i.Classification
]

# clean up working columns
drop_cols <- c(
    "adjacent_var_pos1",
    "adjacent_var_pos2",
    "adjacent_var_pos3",
    "codon_pair_mel",
    "aa_pair_mel",
    "codon_pair_sim",
    "aa_pair_sim",
    "total_codons",
    "shared_codons",
    "Nshared_codons",
    "total_aa",
    "shared_aa",
    "Nshared_aa",
    "mel_aa",
    "sim_aa"
)

shared_dt[, (drop_cols) := NULL]

out_rds <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed.rds")

saveRDS(shared_dt, out_rds)

message("classified RDS written to: ", out_rds)



# CSV SUBSET
subset_table <- shared_dt[1:500, ]
table(subset_table$classification, useNA = "ifany")
fwrite(subset_table, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed_test500.csv")

mkdir -p 
cp /scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed.rds /project/berglandlab/anjali/drosophila_polymorphism/classification

cp /scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/classification/classification_table.csv /project/berglandlab/anjali/drosophila_polymorphism/classification

cp /scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed_test500.csv /project/berglandlab/anjali/drosophila_polymorphism/classification



# ################ CLASSIFICATION - OLD ###############################
# # same-site variants
# shared_dt[
#     ((total_codons == 2 & shared_codons == 2 & Nshared_codons == 0) & 
#     (total_aa == 2 & shared_aa == 2 & Nshared_aa == 0)),
#     classification := "A"
# ]
# # shared_dt[
# #     ((total_codons == 2 & shared_codons == 2 & Nshared_codons == 0) & 
# #     (total_aa == 2 & shared_aa == 2 & Nshared_aa == 0)),
# #     classification := "A"
# # ]

# table(shared_dt$classification, useNA = "ifany")


# ##### OLD: ############################################################

# ##### for same-site (same nt pos) variants, evaluate codon_class on a per-row basis:

# shared_dt[, codon_class := mapply(classify, codon_pair_mel, codon_pair_sim)]
# shared_dt[, aa_class := mapply(classify, aa_pair_mel, aa_pair_sim)]
# table(shared_dt$codon_class, useNA = "ifany")
# table(shared_dt$aa_class, useNA = "ifany")

# # if codon_class = both_shared, then
# #   total_codons = 2, shared_codons = 2, Nshared_codons = 0

# ### add amino acid classification for these classes:
# shared_dt[
#     (codon_class == "both_shared" &
#     aa_class == "both_shared"),
#     classification := "A"
# ]

# # if codon_class = one_shared, then
# #   total_codons = 3, shared_codons = 1, Nshared_codons = 2
# shared_dt[
#     codon_class == "one_shared",
#     classification := "B"
# ]

# # if codon_class = zero_shared, then
# #   total_codons = 4, shared_codons = 0, Nshared_codons = 4
# shared_dt[
#     codon_class == "zero_shared",
#     classification := "C"
# ]

# # shared_dt[, classification := NA_character_]
# table(shared_dt$classification, useNA = "ifany")

# shared_dt[, codon_class := NULL]






# # this input table contains all same-site and different-site (same codon) polymorphisms such that
# #   no more than 1 total polymorphism may exist per codon (only 1 same site OR 1 diff site)

# #### 1. 

# # 2. same site classification: rows where ref_mel and ref_sim are not empty

# # 3. different site classification: 

# # 4. write classification table 
# # same_pos (1 = same, 0 = different) | total_uniqueCodons | shared_codons | Nshared_codons | total_uniqueAA | shared_AA | Nshared_AA

# # save classification table to csv_class
# # fill in classification column in rds, save to rds file




# # ################## testing sets in R
# # a <- unique(sort(c(2,1,3,3)))
# # a
# # b <- unique(sort(c(3,4,1)))
# # b
# # union(a,b) # 1 2 3 4
# # intersect(a,b) # 1 3
# # setdiff(a,b) # 2
# # setdiff(b,a) # 4
# # setequal(a,b) # FALSE