library(data.table)
library(ggplot2)
library(foreach)

# voi <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed_geva.rds")
voi <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed_geva_melMAF5.rds")

load("/scratch/ejy4bu/drosophila/geographic_clines/xtx_c2.Rdata")
xtx <- xc
rm(xc)

load("/scratch/ejy4bu/drosophila/geographic_clines/Drosophila_melanogaster.10_03_2026.nlpTable.Fst.Rdata")
# as nlp

rm(var)
var <- voi[classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y")]
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
    by = classification]

ggplot(data = plot_dt, aes(x=classification, y=mean_XtXst)) + geom_point(size=3) +
    geom_errorbar(aes(ymin = mean_XtXst - stderr, ymax = mean_XtXst + stderr), width = 0.2)

### ANOVA
anova(lm(XtXst ~ class, data=var[!is.na(XtXst)]))
# Analysis of Variance Table

# Response: XtXst
#               Df    Sum Sq Mean Sq F value   Pr(>F)   
# class          1     24805 24804.9  7.2872 0.006946 **
# Residuals 111517 379590629  3403.9                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(lm(XtXst ~ class, data = var[!is.na(XtXst)]))

# Call:
# lm(formula = XtXst ~ class, data = var[!is.na(XtXst)])

# Residuals:
#      Min       1Q   Median       3Q      Max 
# -267.073  -36.366    5.221   42.544  224.840 

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 286.5261     0.3205 894.127  < 2e-16 ***
# classtsp      1.0319     0.3823   2.699  0.00695 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 58.34 on 111517 degrees of freedom
# Multiple R-squared:  6.534e-05,	Adjusted R-squared:  5.638e-05 
# F-statistic: 7.287 on 1 and 111517 DF,  p-value: 0.006946


##### plot XtXst vs geva age estimate 

age <- fread("/scratch/ejy4bu/drosophila/AlleleAges.VA.cm_GEVA.txt")
age[,chr:=tstrsplit(id, "\\.")[[1]]]
age[,pos:=position]

var <- merge(var, age[, .(chr, pos, PostMode, PostMean, PostMedian)], by=c("chr", "pos"), all.x=T)

anova(lm(PostMean ~ class, data = var[!is.na(PostMean)]))

plot_dt <- var[!is.na(PostMean), 
    .(mean_age = mean(PostMean),
    stderr = sd(PostMean)/sqrt(.N),
    n = .N),
    by = classification]

ggplot(data = plot_dt, aes(x=classification, y=mean_age)) + geom_point(size=3) +
    geom_errorbar(aes(ymin = mean_age - stderr, ymax = mean_age + stderr), width = 0.2)

### age vs XtXst 
plot_dt <- var[class=="conv" & !is.na(PostMean) & !is.na(XtXst)]

plot_dt[, age_bin := cut(
    PostMean,
    breaks = seq(min(PostMean), max(PostMean), length.out = 20),
    include.lowest = TRUE
)]

plot_dt <- plot_dt[, .(
    mean_XtXst = mean(XtXst),
    stderr = sd(XtXst)/sqrt(.N),
    n = .N
), by = age_bin]

ggplot(data = plot_dt, aes(x=age_bin, y=mean_XtXst)) + geom_point(size=2) +
    geom_errorbar(aes(ymin = mean_XtXst - stderr, ymax = mean_XtXst + stderr), width = 0.2) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(plot_dt, aes(x = PostMean, y = XtXst)) +
  geom_point(alpha = 0.05, size = 0.3) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_minimal()