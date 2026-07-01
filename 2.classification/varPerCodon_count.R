
library(data.table)

rds_file <- paste0("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_working.rds")
shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")

###############################################################
### Analyzing the Data: Get number polymorphic sites within codons
###############################################################
# this gives a matrix where rows are n_mel meaning 0,1,2,or 3 polymorphic sites within codon
#  and columns are n_sim, similarly

# we are most interested in the class where n_mel and n_sim both = 1, so each species is polymorphic at 1 site within codon
# this will include both same-site (same row) and different site, same codon polymorphisms... 

shared_dt <- unique(shared_dt, by = c("chr", "pos"))

shared_dt[, mel_flag := !is.na(ref_mel)]
shared_dt[, sim_flag := !is.na(ref_sim)]

codon_counts <- shared_dt[,
    .(n_mel = sum(mel_flag), n_sim = sum(sim_flag)),
    by = .(chr, codon_start_pos)
]

codon_counts[n_mel + n_sim > 6]

joint_counts <- codon_counts[, .N, by = .(n_mel, n_sim)]

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

fwrite(joint_matrix,"/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/mel_sim_codonMatrix_NPolymorphicSitesWithinCodon_working.csv")


codon_stats <- shared_dt[
  ,
  .(
    n_mel = sum(!is.na(ref_mel)),
    n_sim = sum(!is.na(ref_sim)),
    n_shared = sum(!is.na(ref_mel) & !is.na(ref_sim))
  ),
  by = .(chr, codon_start_pos)
]

one_one <- codon_stats[n_mel == 1 & n_sim == 1]

one_one[, type := fifelse(n_shared == 1,
                          "same_site",
                          "different_site_within_codon")]

counts <- one_one[, .N, by = type]

message("=== (1,1) codon class breakdown ===")

for (i in seq_len(nrow(counts))) {
  message(counts$type[i], ": ", counts$N[i])
}

# before adding the unique filter:
# n_mel	n_sim_0	n_sim_1	n_sim_2	n_sim_3
# 0	0	0	0	0
# 1	0	185936	1232	6
# 2	0	3165	1331	1
# 3	0	0	0	18

# after adding:
# n_mel	n_sim_0	n_sim_1	n_sim_2	n_sim_3
# 0	0	0	0	0
# 1	0	186274	1163	0
# 2	0	3168	1081	0
# 3	0	0	0	7
