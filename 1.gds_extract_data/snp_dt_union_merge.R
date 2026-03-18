library(SeqArray)
library(data.table)
library(foreach)
library(doMC)

out_dir <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_datatables/test_files/"
out_csv <- paste0(out_dir, "shared_dt_test.csv")
out_rds <- paste0(out_dir, "shared_dt_test.rds")
if(!file.exists(out_csv)) file.create(out_csv)
if(!file.exists(out_rds)) file.create(out_rds)

mel_snp_rds <- paste0(out_dir, "mel_snp_dt.rds")
mel_snp_dt <- readRDS(mel_snp_rds)
sim_snp_rds <- paste0(out_dir, "sim_snp_dt.rds")
sim_snp_dt <- readRDS(sim_snp_dt)

shared <- merge(mel_snp_dt, sim_snp_dt, by = c("chr", "pos", "ref", "alt"), suffixes = c("_mel", "_sim"), all=T)
message(nrow(shared), " total variants")

# shared_test <- shared[1:100]
# shared_table <- get_gds_data(mel_gds, shared, "mel")
message("saving csv to ", out_csv)
fwrite(shared_table, out_csv)
message("saving rds to ", out_rds)
saveRDS(shared_table, out_rds)

message("complete. ", nrow(shared_table), " variants written.")