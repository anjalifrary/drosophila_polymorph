library(SeqArray)
library(data.table)
library(foreach)
library(doMC)

######################################################################

# ### Pool seq files ###

#######################################################################

### INBRED ###
out_dir <- "/scratch/ejy4bu/drosophila/inbred/snpDT/"
in_rds <- paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_clean.rds")
background_rds <- paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.background_1polymorphicCodon.rds")
candidate_rds <- paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.candidate_sharedPolyCodonsOnly.rds")

unfiltered_dt <- readRDS(in_rds)

set.seed(123)

test_dt <- unfiltered_dt[sample(.N, 1000)]  

# nrow(test_dt[!is.na(variant.id_mel)])
# nrow(test_dt[!is.na(variant.id_sim)])
# nrow(test_dt[!is.na(variant.id_mel) & !is.na(variant.id_sim)])