library(data.table)
library(ggplot2)
library(foreach)

voi <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed_geva.rds")

load("/scratch/ejy4bu/drosophila/geographic_clines/xtx_c2.Rdata")
xtx <- xc
rm(xc)

load("/scratch/ejy4bu/drosophila/geographic_clines/Drosophila_melanogaster.10_03_2026.nlpTable.Fst.Rdata")
# as nlp

rm(var)
var <- voi[classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y")]
var <- merge(var, xtx[, .(chr, pos, variant, nLocales_poly, col, global_af, poly_af, XtXst)], by=c("chr", "pos"), all.x=T)
var[classification%in%c("A", "B"), class:="tsp"]
var[classification%in%c("F", "G", "O", "P", "X", "Y"), class:="conv"]

ggplot(data=var, aes(x=group, y=XtXst)) + geom_boxplot()
ggplot(data=var, aes(x=classification, y=XtXst)) + geom_boxplot()
ggplot(data=var, aes(x=group, y=mean(XtXst))) + geom_point()

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
