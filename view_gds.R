library(SeqArray)


### to run:
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1
# Rscript explore_gds.R > /scratch/ejy4bu/drosophila/gds_analysis/view_gds.txt 2>&1

# out file
# output_file <- "/scratch/ejy4bu/drosophila/gds_analysis/view_gds.txt"

# output <- c(
# open files
mel_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.gds")
sim_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds")

cat("\n\n--------Melanogaster------\n\n")
print(seqSummary(mel_gds))
print(ls.gdsn(index.gdsn(mel_gds, path="annotation/info")))

cat("\n\n--------Simulans----------\n\n")
print(seqSummary(sim_gds))
print(ls.gdsn(index.gdsn(sim_gds,path="annotation/info")))

# sample info
print("Mel samples:", head(seqGetData(mel_gds, "sample.id"),10))
print("Sim samples:", head(seqGetData(sim_gds, "sample.id"),10))

# total variants
cat("\n----Number Variants-----\n")
cat("mel variants:", length(seqGetData(mel_gds, "position")), "\n\n")
cat("sim variants:", length(seqGetData(sim_gds, "position")), "\n\n")

# check for annotation field

# check for eff field 

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