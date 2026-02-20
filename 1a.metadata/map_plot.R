# Libraries
library(data.table)
library(maps)
library(ggplot2)


##### A) Map

## set directories
in_dir <- "/project/berglandlab/anjali/drosophila_polymorphism/metadata"
out_dir <- "/project/berglandlab/anjali/drosophila_polymorphism/figures"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

#### samples file
# melanogaster
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
image_file <- file.path(out_dir, "sample_map.png")

png(
  filename = image_file,
  width = 2000,
  height = 1000,
  res = 300
)

print(a)
dev.off()
message(paste0("Saved map to ", image_file))

