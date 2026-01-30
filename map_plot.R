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
samps <- fread(file.path(in_dir,"dest_v2.samps_24Aug2024.xa.csv"))

# simulans
# samps <- fread(file.path(in_dir,"simulans_pooled.meta.use.csv"))

# samps[, set := "All"]
# samps[grepl(":", oldName),set:="Prev. Published"]
# samps[!grepl(":", oldName),set:="Other"]
# formerly Other was named Contaminated DEST - why?

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
    aes(x=long,
        y=lat), 
        size = 2.5, 
        shape = 21
  ) +
  theme_classic(base_size = 14) +
  coord_quickmap() +
  scale_fill_manual(values = c("#f27f65","#8c89c1","#a5c9cc")) +
  scale_color_manual(values = c("#D05438","#615E93","#6F9FA3")) +
  theme(legend.position = "bottom",
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  panel.background = element_blank(), 
  axis.line = element_blank(), 
  axis.text = element_text(color = "black")
  ) 
a

## save file
image_file <- file.path(out_dir, "melanogaster_map.png")

png(
  filename = image_file,
  width = 2000,
  height = 1000,
  res = 300
)

print(a)
dev.off()
message(paste0("Saved map to ", image_file))