# Libraries
library(data.table)
library(maps)
library(ggplot2)


##### A) Map
#### samples file
samps <- fread("/project/berglandlab/anjali/drosophila_polymorphism/metadata/simulans_pooled.meta.use.csv")
samps[grepl(":", oldName),set:="Prev. Published"]
samps[!grepl(":", oldName),set:="Contaminated DEST"]

world <- map_data("world")

a <- ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, map_id = region),
    fill = "#cccccc", 
    color = "darkgray", 
    linewidth = 0.1
  ) + 
  geom_point(
    data = samps,
    aes(x=long,
        y=lat,
        fill = set, color = set), 
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