library(SeqArray)
library(data.table)
library(foreach)
library(doMC)

out_dir <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/"
out_csv <- paste0(out_dir, "all_quality_variants_merge_unfilt_500test.csv")
out_rds <- paste0(out_dir, "all_quality_variants_merge_unfilt.rds")
if(!file.exists(out_csv)) file.create(out_csv)
if(!file.exists(out_rds)) file.create(out_rds)

mel_snp_rds <- paste0(out_dir, "mel_filtered_eff_snp_dt.rds")
mel_snp_dt <- readRDS(mel_snp_rds)
sim_snp_rds <- paste0(out_dir, "sim_filtered_eff_snp_dt.rds")
sim_snp_dt <- readRDS(sim_snp_rds)

# ### test on chromosome 2L
# # filter by 2L chromosome for a smaller subset to test merge on for csv readable
# mel_snp_dt <- mel_snp_dt[chr == "2L"]
# message(nrow(mel_snp_dt), " mel 2L variants")
# sim_snp_dt <- sim_snp_dt[chr == "2L"]
# message(nrow(sim_snp_dt), " sim 2L variants")

# NOTE: this is a union merge... keep ALL variants. no filtering yet
shared_table <- merge(mel_snp_dt, sim_snp_dt, by = c("chr", "pos"), suffixes = c("_mel", "_sim"), all=T)
message(nrow(shared_table), " total variants")

message("saving rds to ", out_rds)
saveRDS(shared_table, out_rds)

message("complete. ", nrow(shared_table), " variants written.")

subset_table <- shared_table[1:500, ]
fwrite(subset_table, out_csv)
message("saved first 500 rows to csv at ", out_csv)
