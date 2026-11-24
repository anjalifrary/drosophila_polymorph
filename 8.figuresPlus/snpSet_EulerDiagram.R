library(data.table)

# load tables (including non-shared and shared variation, both TSP/conv and different)
poolseq_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/poolseq_dest/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.shared.bothMelSim.classed.rds")
inbred_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.11_2026.shared.bothMelSim.0.05MAFfilter.classed.rds")


# Classification groups
tsp_classes <- c("A", "B")
conv_classes <- c("F", "G", "O", "P", "X", "Y")

get_euler_counts <- function(dt, dataset_name) {

  # Species candidate sets
  mel <- unique(dt[!is.na(variant.id_mel), .(chr, pos)])
  sim <- unique(dt[!is.na(variant.id_sim), .(chr, pos)])

  # Shared variation
  shared <- dt[!is.na(variant.id_mel) & !is.na(variant.id_sim)]

  # Classification groups
  shared[classification %in% tsp_classes, group := "TSP"]
  shared[classification %in% conv_classes, group := "Convergent"]
  shared[!is.na(classification) & !classification %in%c(tsp_classes, conv_classes), group := "Shared other"]

  data.table(
    dataset = dataset_name,
    mel = nrow(mel),
    sim = nrow(sim),
    shared = nrow(shared),
    TSP = nrow(shared[group == "TSP"]),
    convergent = nrow(shared[group == "Convergent"]),
    other_shared = nrow(shared[group == "Shared Other"]),
    mel_only = nrow(mel) - nrow(shared), 
    sim_only = nrow(sim) - nrow(shared)
  )
}

euler_counts <- rbind(
  get_euler_counts(poolseq_dt, "PoolSeq"),
  get_euler_counts(inbred_dt, "Inbred")
)



#######################################################################
### Make Diagrams

## had to load this in local Rstudio 
library(eulerr)
library(data.table)

euler_counts <- fread("~/UVA/Bergland/euler_counts.csv")

make_euler <- function(dataset_name) {
    d <- euler_counts[dataset == dataset_name]

    eu <- euler(c(
        "D. simulans" = d$sim_only,
        "D. melanogaster" = d$mel_only,
        "D. simulans&D. melanogaster" = d$shared
    ))

    if (dataset_name == "Inbred") {
        g <- plot(
            eu,
            quantities = FALSE,
            labels = FALSE,
            fills = c("#9ad9db", "#f49892"),
            edges = list(col = "grey20", lwd = 1),
            main = "Inbred"
        )
        print(g)

        # Species-only counts
        grid::grid.text(prettyNum(d$sim_only, big.mark = ","),
            x = 0.30, y = 0.06, gp = grid::gpar(fontsize = 12))
        grid::grid.text(prettyNum(d$mel_only, big.mark = ","),
            x = 0.70, y = 0.06, gp = grid::gpar(fontsize = 12))

        # Intersection count
        grid::grid.lines(x = c(0.572, 0.572), y = c(0.43, 0.25), gp = grid::gpar(col = "grey20", lwd = 1, lty = "dotted"))
        grid::grid.text(prettyNum(d$shared, big.mark = ","), x = 0.572, y = 0.19, gp = grid::gpar(fontsize = 12))

        # Species labels
        grid::grid.text("D. simulans", x = 0.25, y = 0.12, gp = grid::gpar(fontsize = 12))
        grid::grid.text("D. melanogaster", x = 0.75, y = 0.12, gp = grid::gpar(fontsize = 12))

    } else if (dataset_name == "PoolSeq") {

        g <- plot(
            eu,
            quantities = FALSE,
            labels = FALSE,
            fills = c("#9ad9db", "#f49892"),
            edges = list(col = "grey20", lwd = 1),
            main = "PoolSeq"
        )
        print(g)

        # number labels 
        grid::grid.text(
            prettyNum(d$sim_only, big.mark = ","),
            x = 0.25, y = 0.06,
            gp = grid::gpar(fontsize = 12)
        )

        grid::grid.text(
            prettyNum(d$mel_only, big.mark = ","),
            x = 0.75, y = 0.06,
            gp = grid::gpar(fontsize = 12)
        )

        # Interesection count + dotted line
        grid::grid.lines(
            x = c(0.445, 0.445),
            y = c(0.43, 0.25),
            gp = grid::gpar(
                col = "grey20",
                lwd = 1,
                lty = "dotted"
            )
        )
        grid::grid.text(
            prettyNum(d$shared, big.mark = ","),
            x = 0.445, y = 0.19,
            gp = grid::gpar(fontsize = 12)
        )
        
        # Species labels beneath circles
        grid::grid.text(
            "D. simulans",
            x = 0.25, y = 0.12,
            gp = grid::gpar(fontsize = 12)
        )

        grid::grid.text(
            "D. melanogaster",
            x = 0.75, y = 0.12,
            gp = grid::gpar(fontsize = 12)
        )
    }
}

make_euler("PoolSeq")
make_euler("Inbred")

plot_composition <- function(dataset_name) {

    ggplot(
        shared_composition_long[dataset == dataset_name],
        aes(x = percent, y = dataset, fill = class)
    ) +
        geom_col(width = 0.35) +
        geom_text(
            aes(label = ifelse(
                percent >= 0.08,
                paste0(round(percent * 100, 1), "%"),
                ""
            )),
            position = position_stack(vjust = 0.5),
            size = 2.5
        ) +
        scale_fill_manual(
            values = c(
                "TSP" = "#8fb3dc",
                "Conv" = "#9ed6a6",
                "Unclassified" = "#D9D9D9"
            ),
            labels = c("TSP candidates", "Convergent", "Unclassified")
        ) +
        scale_x_continuous(
            labels = scales::percent,
            expand = c(0, 0)
        ) +
        labs(
            x = NULL,
            y = NULL,
            fill = NULL
        ) +
        theme_void() +
        theme(
            legend.position = "bottom",
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.25, "cm"),
            legend.spacing.x = unit(0.1, "cm"),
            plot.margin = margin(0, 2, 0, 2)
        )
}

plot_composition("PoolSeq")
plot_composition("Inbred")




library(ggplot2)
library(ggforce)

### too cartoon & colors don't match
plot_euler <- function(counts, dataset_name) {

  d <- counts[dataset == dataset_name]

  sim_only <- d$sim - d$shared
  mel_only <- d$mel - d$shared

  ggplot() +

    # Species circles
    geom_circle(
      aes(x0 = -1, y0 = 0, r = 2),
      fill = "#F8766D", alpha = 0.20,
      color = "#F8766D", linewidth = 0.8
    ) +
    geom_circle(
      aes(x0 = 1, y0 = 0, r = 2),
      fill = "#00BFC4", alpha = 0.20,
      color = "#00BFC4", linewidth = 0.8
    ) +

    # Species labels
    annotate("text", x = -1.5, y = 1.55,
             label = "D. simulans", fontface = "italic", size = 5) +
    annotate("text", x = 1.5, y = 1.55,
             label = "D. melanogaster", fontface = "italic", size = 5) +

    # Species-only counts
    annotate("text", x = -1.8, y = 0,
             label = format(sim_only, big.mark = ","),
             size = 5) +
    annotate("text", x = 1.8, y = 0,
             label = format(mel_only, big.mark = ","),
             size = 5) +

    # Shared classification breakdown
    annotate("text", x = 0, y = 0.65,
             label = paste0("Shared\n", format(d$shared, big.mark = ",")),
             fontface = "bold", size = 4.5) +
    annotate("text", x = 0, y = 0.05,
             label = paste0("TSP: ", format(d$tsp, big.mark = ",")),
             size = 4) +
    annotate("text", x = 0, y = -0.4,
             label = paste0("Convergent: ", format(d$convergent, big.mark = ",")),
             size = 3.8) +
    annotate("text", x = 0, y = -0.8,
             label = paste0("Other: ", format(d$other_shared, big.mark = ",")),
             size = 3.8) +

    labs(title = dataset_name) +

    coord_fixed(xlim = c(-3.2, 3.2), ylim = c(-2.5, 2.5)) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.margin = margin(15, 15, 15, 15)
    )
}
plot_euler(euler_counts, "PoolSeq")
plot_euler(euler_counts, "Inbred")


### also ugly:
library(VennDiagram)
library(grid)

draw.pairwise.venn(
  area1 = euler_counts[dataset == "PoolSeq", sim],
  area2 = euler_counts[dataset == "PoolSeq", mel],
  cross.area = euler_counts[dataset == "PoolSeq", shared],
  category = c("D. simulans", "D. melanogaster"),
  main = "PoolSeq",
  scaled = TRUE
)

draw.pairwise.venn(
  area1 = euler_counts[dataset == "Inbred", sim],
  area2 = euler_counts[dataset == "Inbred", mel],
  cross.area = euler_counts[dataset == "Inbred", shared],
  category = c("D. simulans", "D. melanogaster"),
  main = "Inbred",
  scaled = TRUE
)


plot_euler(euler_counts, "PoolSeq")

plot_euler(euler_counts, "Inbred")










### UGLY
library(ggplot2)

# Function to draw circles
circle <- function(x, y, r) {
  angle <- seq(0, 2*pi, length.out = 200)
  data.frame(x = x + r*cos(angle), y = y + r*sin(angle))
}

plot_euler <- function(counts, dataset_name) {

  d <- counts[dataset == dataset_name]

  # Outer circles
  sim_circle <- circle(-0.85, 0, 1.6)
  mel_circle <- circle( 0.85, 0, 1.6)

  # Smaller circles inside shared overlap
  tsp_circle  <- circle(0, 0.35, 0.22)
  conv_circle <- circle(0, -0.35, 0.22)

  ggplot() +

    # Species circles: default ggplot colors
    geom_path(data = sim_circle, aes(x, y, color = "D. simulans"),
              linewidth = 1) +
    geom_path(data = mel_circle, aes(x, y, color = "D. melanogaster"),
              linewidth = 1) +

    # Subset circles
    geom_path(data = tsp_circle, aes(x, y, color = "TSP"),
              linewidth = 0.8) +
    geom_path(data = conv_circle, aes(x, y, color = "Convergent"),
              linewidth = 0.8) +

    # Species counts
    annotate("text", x = -1.3, y = 1.85,
             label = paste0("D. simulans\n", format(d$sim, big.mark = ",")),
             size = 4) +
    annotate("text", x = 1.3, y = 1.85,
             label = paste0("D. melanogaster\n", format(d$mel, big.mark = ",")),
             size = 4) +

    # Shared count
    annotate("text", x = 0, y = 0.95,
             label = paste0("Shared\n", format(d$shared, big.mark = ",")),
             size = 4) +

    # TSP and convergent counts
    annotate("text", x = 0, y = 0.35,
             label = paste0("TSP\n", format(d$tsp, big.mark = ",")),
             size = 3) +
    annotate("text", x = 0, y = -0.35,
             label = paste0("Convergent\n", format(d$convergent, big.mark = ",")),
             size = 2.8) +

    scale_color_discrete() +

    coord_fixed(xlim = c(-2.8, 2.8), ylim = c(-2.2, 2.3)) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    ) +
    labs(title = dataset_name)
}

plot_euler(euler_counts, "PoolSeq")

plot_euler(euler_counts, "Inbred")