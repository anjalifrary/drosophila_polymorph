# Libraries
library(data.table)
library(maps)
library(ggplot2)
library(SeqArray)

##### A) Map

## set directories
in_dir <- "/project/berglandlab/anjali/drosophila_polymorphism/metadata"
out_dir <- "/scratch/ejy4bu/drosophila/gds_analysis/figures/map"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

## Load metadata files
mel_metadata <- fread(file.path(in_dir,"dest_v2.samps_24Aug2024.xa.csv"))
mel_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.eff.gds")
sim_metadata <- fread(file.path(in_dir,"simulans_pooled.meta.use.csv"))
sim_contam <- fread(file.path(in_dir, "Dsim_contamRates.csv"))
sim_gds <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.eff.gds")


# # Get samples that pass quality filter 

# mel_metadata[, species :="D. melanogaster"]
# good_mel_samps <- mel_metadata[Recommendation == "Pass", sampleId]
# good_mel_samps <- intersect(seqGetData(mel_gds, "sample.id"), good_mel_samps) 
# # idk why the gds doesn't contain these samples but 530 samples without merging and 502 samples with merging
# message(length(good_mel_samps), " = 502 mel samples pass filter")

# sim_metadata[, species :="D. simulans"]
# good_sim_samps <- sim_contam[propSim >= 0.95, sampleId]
# good_sim_samps <- intersect(seqGetData(sim_gds, "sample.id"), good_sim_samps)
# message(length(good_sim_samps), " = 15 sim samples pass filter")

# #combine data sets
# mel_dt <- mel_metadata[sampleId %in% good_mel_samps]
# sim_dt <- sim_metadata[identifier %in% good_sim_samps]
# samps <- rbindlist(list(mel_dt, sim_dt), fill = TRUE)

## For unfiltered map:
mel <- fread(file.path(in_dir,"dest_v2.samps_24Aug2024.xa.csv"))
mel[, species :="D. melanogaster"]

# simulans
sim <- fread(file.path(in_dir,"simulans_pooled.meta.use.csv"))
sim[, species :="D. simulans"]

#combine data sets
samps <- rbindlist(list(mel, sim), fill = TRUE)

world <- map_data("world")

a <- ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group=group),
    fill = "#cccccc", 
    color = "darkgray", 
    linewidth = 0.1
  ) + 
  geom_point(
    data = samps,
    aes(x=long, y=lat, fill = species), 
    color = "black",
        size = 2.5, 
        shape = 21
        # fill = "#cccccc", color = "#cccccc"
  ) +
  theme_classic(base_size = 12) +
  coord_quickmap() +
  scale_fill_manual(values = c("D. melanogaster"="#f27f65", "D. simulans"="#8c89c1")) +
  # scale_color_manual(values = c("D. mel"="#D05438", "D. sim"="#615E93")) +
  theme(legend.position = "bottom",
  panel.grid = element_blank(), 
  # panel.background = element_blank(), 
  axis.line = element_blank(), 
  axis.text = element_text(color = "black", size=8)
  ) 
a

### Save file
image_file <- file.path(out_dir, "sample_unfiltered_map.png")

png(
  filename = image_file,
  width = 2000,
  height = 1000,
  res = 300
)

print(a)
dev.off()
message(paste0("Saved map to ", image_file))
