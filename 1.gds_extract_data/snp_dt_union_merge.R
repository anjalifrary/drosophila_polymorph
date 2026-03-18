

shared <- merge(mel_snp_dt, sim_snp_dt, by = c("chr", "pos"), suffixes = c("_mel", "_sim"), all=T)
message(nrow(shared), " total variants")

# shared_test <- shared[1:100]
shared_table <- get_gds_data(mel_gds, shared, "mel")
# message("saving csv to ", out_csv)
# fwrite(shared_table, out_csv)
message("saving rds to ", out_rds)
saveRDS(shared_table, out_rds)

message("complete. ", nrow(shared_table), " variants written.")