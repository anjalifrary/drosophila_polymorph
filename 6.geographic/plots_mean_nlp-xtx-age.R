library(data.table)
library(ggplot2)
library(foreach)

# voi <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed_geva.rds")
voi <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_MAF5_06-18-2026.rds")

load("/scratch/ejy4bu/drosophila/geographic_clines/xtx_c2.Rdata")
xtx <- xc
rm(xc)

load("/scratch/ejy4bu/drosophila/geographic_clines/Drosophila_melanogaster.10_03_2026.nlpTable.Fst.Rdata")
mel_nlp <- nlp  
rm(nlp)
load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_simulans.17_06_2026.nlpTable.Rdata")
sim_nlp <- nlp
rm(nlp)

rm(var)
cand_classes <- c("A", "B", "F", "G", "O", "P", "X", "Y")
var <- voi[classification%in%(cand_classes)]
var <- merge(var, xtx[, .(chr, pos, col, XtXst)], by=c("chr", "pos"), all.x=T)
var[classification%in%c("A", "B"), class:="tsp"]
var[classification%in%c("F", "G", "O", "P", "X", "Y"), class:="conv"]

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

age <- fread("/scratch/ejy4bu/drosophila/AlleleAges.VA.cm_GEVA.txt")
age[,chr:=tstrsplit(id, "\\.")[[1]]]
age[,pos:=position]

var <- merge(var, age[, .(chr, pos, PostMode, PostMean, PostMedian)], by=c("chr", "pos"), all.x=T)

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
    by = classification]

ggplot(data = plot_dt, aes(x=classification, y=mean_nlp)) + geom_point(size=3) +
    geom_errorbar(aes(ymin = mean_nlp - stderr, ymax = mean_nlp + stderr), width = 0.2)


# plot mel and sim mean nlp color-coded 
plot_dt <- rbindlist(list(
  var[!is.na(nLocales_poly_mel),
      .(mean_nlp = mean(nLocales_poly_mel),
        stderr = sd(nLocales_poly_mel)/sqrt(.N),
        n = .N,
        species = "mel"),
      by = class],

  var[!is.na(nLocales_poly_sim),
      .(mean_nlp = mean(nLocales_poly_sim),
        stderr = sd(nLocales_poly_sim)/sqrt(.N),
        n = .N,
        species = "sim"),
      by = class]
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

var_filt <- var[nLocales_poly_mel>10 & nLocales_poly_sim>10]

plot_dt <- rbindlist(lapply(p_vals, function(p) {

  tmp <- get_tail_props(var_filt, p)

  melt(tmp,
       id.vars = "class",
       measure.vars = c("lower_tail", "upper_tail"),
       variable.name = "tail",
       value.name = "prop")[, `:=`(p = p)]
}))

ggplot(plot_dt, aes(x = class, y = prop, fill = tail)) +
  geom_col(position = "dodge") +
  facet_wrap(~p, labeller = label_both)


