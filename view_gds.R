library(SeqArray)

mel_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.gds")
sim_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds")

# Just print everything and redirect with > in bash
print(seqSummary(mel_gds))
print(seqSummary(sim_gds))
print(ls.gdsn(index.gdsn(mel_gds, "annotation/info")))
print(ls.gdsn(index.gdsn(sim_gds, "annotation/info")))

seqClose(mel_gds)
seqClose(sim_gds)

# # from claude

# # Try to get annotation data
# try({
#   ann_test <- seqGetData(mel_gds, "annotation/info/ANN")
#   print("Found ANN field!")
#   print(head(ann_test))
# }, silent = FALSE)

# # If no ANN, check for other annotation fields
# try({
#   eff_test <- seqGetData(mel_gds, "annotation/info/EFF")
#   print("Found EFF field!")
#   print(head(eff_test))
# }, silent = FALSE)