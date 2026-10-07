library(data.table)
library(ggplot2)
library(foreach)
### old 
# voi <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed_geva.rds")
# voi <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_MAF5_06-18-2026.rds")

### inbred:
# voi <- readRDS("/scratch/ejy4bu/drosophila/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.candidatesABFGOPXY.classed.rds")
inbred <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.candidatesABFGOPXY.classed.MAF.rds")

### poolseq:
# voi <- readRDS("/scratch/ejy4bu/drosophila/DEST_remake/snpDT/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.candidatesABFGOPXY.classed.rds")
poolseq <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/poolseq_dest/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.candidatesABFGOPXY.classed.rds")

load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/xtx_c2.Rdata")
xtx <- xc
rm(xc)

load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_melanogaster.11_08_2026.nlpTable.paramask.genmap.busco.repeatmask.wmdust.Rdata")
mel_nlp <- nlp  
rm(nlp)
load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_simulans.17_06_2026.nlpTable.Rdata")
sim_nlp <- nlp
rm(nlp)


age <- fread("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/AlleleAges.VA.cm_GEVA.txt")
age[,chr:=tstrsplit(id, "\\.")[[1]]]
age[,pos:=position]

cand_classes <- c("A", "B", "F", "G", "O", "P", "X", "Y")


# rm(var)
cand_classes <- c("A", "B", "F", "G", "O", "P", "X", "Y")
var <- voi[classification%in%(cand_classes)]
var <- merge(var, xtx[, .(chr, pos, XtXst)], by.y=c("chr", "pos"), by.x=c("chr_dm6", "pos_dm6"), all.x=T)
var[classification%in%c("A", "B"), class:="tsp"]
var[classification%in%c("F", "G", "O", "P", "X", "Y"), class:="conv"]

var <- merge(var, mel_nlp[, .(chr, pos, nLocales_poly)],
    by.x=c("chr_dm6", "pos_dm6"), by.y=c("chr", "pos"), all.x=T)
setnames(var, "nLocales_poly", "nLocales_poly_mel")

var <- merge(var, sim_nlp[, .(chr, pos, nLocales_poly)], 
    by.x=c("chr_dm6", "pos_dm6"), by.y=c("chr", "pos"),  all.x=T)
setnames(var, "nLocales_poly", "nLocales_poly_sim")



var <- merge(var, age[, .(chr, pos, PostMode, PostMean, PostMedian)], by.y=c("chr", "pos"), by.x=c("chr_dm6", "pos_dm6"), all.x=T)


ggplot(data=var, aes(x=class, y=XtXst)) + geom_boxplot()
ggplot(data=var, aes(x=classification, y=XtXst)) + geom_boxplot()
# ggplot(data=var, aes(x=class, y=mean(XtXst))) + geom_point()

plot_dt <- var[!is.na(XtXst), 
    .(mean_XtXst = mean(XtXst),
    stderr = sd(XtXst)/sqrt(.N),
    n = .N),
    by = class]

ggplot(data = plot_dt, aes(x=class, y=mean_XtXst)) + geom_point(size=3) +
    geom_errorbar(aes(ymin = mean_XtXst - stderr, ymax = mean_XtXst + stderr), width = 0.2)

### ANOVA
anova(lm(XtXst ~ class, data=var[!is.na(XtXst)]))
summary(lm(XtXst ~ class, data = var[!is.na(XtXst)]))

### age from geva - mel-only age estimates 

anova(lm(PostMean ~ class, data = var[!is.na(PostMean)]))

plot_dt <- var[!is.na(PostMean), 
    .(mean_age = mean(PostMean),
    stderr = sd(PostMean)/sqrt(.N),
    n = .N),
    by = class]

ggplot(data = plot_dt, aes(x=class, y=mean_age)) + geom_point(size=3) +
    geom_errorbar(aes(ymin = mean_age - stderr, ymax = mean_age + stderr), width = 0.2)

# ### skip plot, too confusing
# ### age vs XtXst 
# plot_dt <- var[class=="conv" & !is.na(PostMean) & !is.na(XtXst)]

# plot_dt[, age_bin := cut(
#     PostMean,
#     breaks = seq(min(PostMean), max(PostMean), length.out = 20),
#     include.lowest = TRUE
# )]

# plot_dt <- plot_dt[, .(
#     mean_XtXst = mean(XtXst),
#     stderr = sd(XtXst)/sqrt(.N),
#     n = .N
# ), by = age_bin]

# ggplot(data = plot_dt, aes(x=age_bin, y=mean_XtXst)) + geom_point(size=2) +
#     geom_errorbar(aes(ymin = mean_XtXst - stderr, ymax = mean_XtXst + stderr), width = 0.2) + 
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ggplot(plot_dt, aes(x = PostMean, y = XtXst)) +
#   geom_point(alpha = 0.05, size = 0.3) +
#   geom_smooth(method = "loess", se = TRUE) +
#   theme_minimal()



########################################################3
# mean nlp vs class(ification)

spp = "mel"
# spp = "sim"

nlp <- paste0("nLocales_poly_", spp)

plot_dt <- var[!is.na(get(nlp)), 
    .(mean_nlp = mean(get(nlp)),
    stderr = sd(get(nlp))/sqrt(.N),
    n = .N),
    by = class]

ggplot(data = plot_dt, aes(x=class, y=mean_nlp)) + geom_point(size=3) +
    geom_errorbar(aes(ymin = mean_nlp - stderr, ymax = mean_nlp + stderr), width = 0.2)


# plot mel and sim mean nlp color-coded 
plot_dt <- rbindlist(list(
  var[!is.na(nLocales_poly_mel),
      .(mean_nlp = mean(nLocales_poly_mel),
        stderr = sd(nLocales_poly_mel)/sqrt(.N),
        n = .N,
        species = "mel"),
      by = class]
#       ,

#   var[!is.na(nLocales_poly_sim),
#       .(mean_nlp = mean(nLocales_poly_sim),
#         stderr = sd(nLocales_poly_sim)/sqrt(.N),
#         n = .N,
#         species = "sim"),
#       by = class]
))

ggplot(plot_dt, aes(x = class, y = mean_nlp, color = species, group = species)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = mean_nlp - stderr, ymax = mean_nlp + stderr),
    width = 0.2, position = position_dodge(width = 0.5))


################################################3
### XtXst - tails
# estimate proportion of conv/tsp in tails
# filter for MAF, global, or nlp frequency 

# what fraction of conv variants are in the lower vs upper tail? tsp variants? 
    # can't ask "in the tail, what fraction are conv vs tsp" because num tsp >> conv

get_tail_props <- function(data, p = 0.05) {

  lower <- quantile(data$XtXst, p, na.rm = TRUE)
  upper <- quantile(data$XtXst, 1 - p, na.rm = TRUE)

  data[, .(
    lower_tail  = mean(XtXst <= lower, na.rm = TRUE),
    upper_tail = mean(XtXst >= upper, na.rm = TRUE)
  ), by = class]
}

p_vals <- c(
    0.01,
    0.05,
    0.10,
    0.30
)

var_filt <- var[
    nLocales_poly_mel > 20 & !is.na(XtXst)
    # nLocales_poly_sim > 20
]

plot_dt <- rbindlist(lapply(p_vals, function(p) {

  tmp <- get_tail_props(var_filt, p)

  melt(tmp,
       id.vars = "class",
       measure.vars = c("lower_tail", "upper_tail"),
       variable.name = "tail",
       value.name = "prop")[, `:=`
       (p = p,
       prop_norm = prop / p
       )]
}))

ggplot(plot_dt, aes(x = class, y = prop_norm, fill = tail)) +
  geom_col(position = "dodge") +
  facet_wrap(~p, labeller = label_both) +
  geom_hline(
        yintercept = 1,
        linetype = "dashed"
    )




### XtXst: proportion of TSP / Conv variants in tail as a function of p
### mel NLP only

var_filt <- var[
  nLocales_poly_mel > 20 &
  !is.na(XtXst)
]

get_tail_props <- function(data, p) {

  lower <- quantile(data$XtXst, p, na.rm = TRUE)
  upper <- quantile(data$XtXst, 1 - p, na.rm = TRUE)

  data[, .(
    lower_tail = mean(XtXst <= lower),
    upper_tail = mean(XtXst >= upper)
  ), by = class]
}

p_vals <- seq(0.01, 0.50, by = 0.01)

plot_dt <- rbindlist(lapply(p_vals, function(p) {

  tmp <- get_tail_props(var_filt, p)

  melt(
    tmp,
    id.vars = "class",
    measure.vars = c("lower_tail", "upper_tail"),
    variable.name = "tail",
    value.name = "prop"
  )[, p := p]
}))

plot_dt[, tail := fifelse(
  tail == "lower_tail",
  "Lower",
  "Upper"
)]

ggplot(
  plot_dt,
  aes(
    x = p,
    y = prop,
    color = class,
    linetype = tail,
    group = interaction(class, tail)
  )
) +
  geom_line(linewidth = 1) +
  geom_hline(
    aes(yintercept = p),
    linetype = "dashed",
    color = "grey50"
  ) +
  scale_x_continuous(
    labels = scales::percent,
    breaks = c(.01, .05, .10, .25, .50)
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    x = "Tail size (p)",
    y = "Proportion of variants in tail",
    color = NULL,
    linetype = NULL
  ) +
  theme_classic()

  ### y-axis = enrichment
plot_dt[, enrichment := prop / p]

ggplot(
  plot_dt,
  aes(
    x = p,
    y = enrichment,
    color = class,
    linetype = tail,
    group = interaction(class, tail)
  )
) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_x_continuous(
    labels = scales::percent,
    breaks = c(.01, .05, .10, .25, .50)
  ) +
  labs(
    x = "Tail size (p)",
    y = "Tail enrichment",
    color = NULL,
    linetype = NULL
  ) +
  theme_classic()




### INBRED & POOLSEQ
prepare_var <- function(voi) {

  var <- voi[classification %in% cand_classes]

  # XtXst
  var <- merge(
    var,
    xtx[, .(chr, pos, XtXst)],
    by.y = c("chr", "pos"),
    by.x = c("chr_dm6", "pos_dm6"),
    all.x = TRUE
  )

  # TSP vs Conv
  var[
    classification %in% c("A", "B"),
    class := "TSP"
  ]

  var[
    classification %in% c("F", "G", "O", "P", "X", "Y"),
    class := "Conv"
  ]

  # mel NLP ONLY
  var <- merge(
    var,
    mel_nlp[, .(chr, pos, nLocales_poly)],
    by.x = c("chr_dm6", "pos_dm6"),
    by.y = c("chr", "pos"),
    all.x = TRUE
  )

  setnames(var, "nLocales_poly", "nLocales_poly_mel")

  # GEVA age
  var <- merge(
    var,
    age[, .(chr, pos, PostMode, PostMean, PostMedian)],
    by.y = c("chr", "pos"),
    by.x = c("chr_dm6", "pos_dm6"),
    all.x = TRUE
  )

  var
}

inbred_var  <- prepare_var(inbred)
poolseq_var <- prepare_var(poolseq)

### AGE:
age_dt <- rbind(
  inbred_var[, .(
    mean_age = mean(PostMean, na.rm = TRUE),
    stderr = sd(PostMean, na.rm = TRUE) / sqrt(sum(!is.na(PostMean)))
  ), by = class][, dataset := "Inbred"],

  poolseq_var[, .(
    mean_age = mean(PostMean, na.rm = TRUE),
    stderr = sd(PostMean, na.rm = TRUE) / sqrt(sum(!is.na(PostMean)))
  ), by = class][, dataset := "PoolSeq"]
)
p_age <- ggplot(
  age_dt,
  aes(x = class, y = mean_age, color = dataset)
) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 3
  ) +
  geom_errorbar(
    aes(
      ymin = mean_age - stderr,
      ymax = mean_age + stderr
    ),
    position = position_dodge(width = 0.4),
    width = 0.2
  ) +
  labs(
    x = NULL,
    y = "Mean allele age",
    color = NULL
  ) +
  theme_classic()

### NLP
nlp_dt <- rbind(
  inbred_var[nLocales_poly_mel > 20, .(
    mean_nlp = mean(nLocales_poly_mel, na.rm = TRUE),
    stderr = sd(nLocales_poly_mel, na.rm = TRUE) /
      sqrt(sum(!is.na(nLocales_poly_mel)))
  ), by = class][, dataset := "Inbred"],

  poolseq_var[nLocales_poly_mel > 20, .(
    mean_nlp = mean(nLocales_poly_mel, na.rm = TRUE),
    stderr = sd(nLocales_poly_mel, na.rm = TRUE) /
      sqrt(sum(!is.na(nLocales_poly_mel)))
  ), by = class][, dataset := "PoolSeq"]
)
p_nlp <- ggplot(
  nlp_dt,
  aes(x = class, y = mean_nlp, color = dataset)
) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 3
  ) +
  geom_errorbar(
    aes(
      ymin = mean_nlp - stderr,
      ymax = mean_nlp + stderr
    ),
    position = position_dodge(width = 0.4),
    width = 0.2
  ) +
  labs(
    x = NULL,
    y = "Mean mel NLP",
    color = NULL
  ) +
  theme_classic()

### XTX
get_tail_enrichment <- function(data, dataset) {

  data <- data[
    nLocales_poly_mel > 20 &
    !is.na(XtXst)
  ]

  p_vals <- seq(.01, .50, .01)

  rbindlist(lapply(p_vals, function(p) {

    lower <- quantile(data$XtXst, p, na.rm = TRUE)
    upper <- quantile(data$XtXst, 1 - p, na.rm = TRUE)

    rbind(
      data.table(
        p = p,
        class = "TSP",
        tail = "Lower",
        enrichment = mean(data[class == "TSP", XtXst] <= lower) / p
      ),
      data.table(
        p = p,
        class = "TSP",
        tail = "Upper",
        enrichment = mean(data[class == "TSP", XtXst] >= upper) / p
      ),
      data.table(
        p = p,
        class = "Conv",
        tail = "Lower",
        enrichment = mean(data[class == "Conv", XtXst] <= lower) / p
      ),
      data.table(
        p = p,
        class = "Conv",
        tail = "Upper",
        enrichment = mean(data[class == "Conv", XtXst] >= upper) / p
      )
    )

  }))[, dataset := dataset]
}


tail_dt <- rbind(
  get_tail_enrichment(inbred_var, "Inbred"),
  get_tail_enrichment(poolseq_var, "PoolSeq")
)

p_xtx_tail <- ggplot(
  tail_dt,
  aes(
    x = p,
    y = enrichment,
    color = dataset,
    linetype = class,
    group = interaction(dataset, class)
  )
) +
  geom_line(linewidth = 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  scale_linetype_manual(
    values = c(
      "TSP" = "solid",
      "Conv" = "11"
    )) +
  facet_wrap(~tail) +
  scale_x_continuous(
    labels = scales::percent,
    breaks = c(.01, .05, .10, .25, .50)
  ) +
  labs(
    x = "Tail size",
    y = "XtXst tail enrichment",
    color = NULL,
    linetype = NULL
  ) +
  theme_classic()

# smoother
p_xtx_tail <- ggplot(
  tail_dt,
  aes(
    x = p,
    y = enrichment,
    color = dataset,
    linetype = class,
    group = interaction(dataset, class)
  )
) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    span = 0.35,
    linewidth = 0.9
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey60"
  ) +
  facet_wrap(~tail, nrow = 1) +
  scale_x_continuous(
    labels = scales::percent,
    breaks = c(.01, .05, .10, .25, .50)
  ) +
  scale_linetype_manual(
    values = c(
      "TSP" = "solid",
      "Conv" = "dashed"
    )
  ) +
  labs(
    x = "Tail size",
    y = "XtXst tail enrichment",
    color = NULL,
    linetype = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# move legend 
p_xtx_tail <- ggplot(
  tail_dt,
  aes(
    x = p,
    y = enrichment,
    color = dataset,
    linetype = class,
    group = interaction(dataset, class)
  )
) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    span = 0.35,
    linewidth = 0.9
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey60"
  ) +
  facet_wrap(~tail, nrow = 1) +
  scale_x_continuous(
  limits = c(.01, .50),
  labels = scales::percent,
  breaks = c(.10, .25, .50)
) + 
  scale_linetype_manual(
    values = c(
      "TSP" = "solid",
      "Conv" = "33"
    )
  ) +
  labs(
    x = "Tail size",
    y = "XtXst tail enrichment",
    color = NULL,
    linetype = NULL
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(
      order = 2,
      override.aes = list(color = "grey20")
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

library(patchwork)
### 3 PANELS
# side by side by side
figure1 <- p_age + p_xtx_tail + p_nlp +
  plot_layout(widths = c(1, 3, 1))

# age and nlp side by side top
# xtx wide on bottom
p_age2 <- p_age +
  theme(legend.position = "none")

p_nlp2 <- p_nlp +
  theme(legend.position = "none")

p_xtx_tail2 <- p_xtx_tail +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(
      order = 2,
      override.aes = list(color = "grey20")
    )
  ) +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )

  figure1 <- (
  p_age2 + p_nlp2
) / p_xtx_tail2 +
  plot_layout(
    heights = c(1, 1.25)
  ) +
  plot_annotation(tag_levels = "A")