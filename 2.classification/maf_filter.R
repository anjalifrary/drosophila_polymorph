library(data.table)
library(ggplot2)
library(foreach)
library(SeqArray)

# sim file
load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_simulans.17_06_2026.nlpTable.Rdata")
# load("/project/berglandlab/multispecies_endemism/data/collectiveAnalysis_version3/Drosophila_simulans.10_03_2026.nlpTable.paramask.genmap.busco.Rdata")
### this sim file is in simulans genospace, wrong coordinates
sim_nlp <- nlp 
rm(nlp)

# mel file
# load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_melanogaster.17_06_2026.nlpTable.Rdata")
load("/project/berglandlab/multispecies_endemism/data/collectiveAnalysis_version3/Drosophila_melanogaster.10_03_2026.nlpTable.paramask.genmap.busco.Rdata")
# confirm these nlp files are the same... (excluding extra columns)
mel_nlp <- nlp
rm(nlp)   

# all variant table, no maf filter yet:
variants <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/noMAFfilter/subset_qualVar_ofInterest.rds")
sim_nlp <- merge(sim_nlp, variants[, c("chr", "pos", "codon_start_pos", "classification")], by=c("chr", "pos"), all.x=T)
mel_nlp <- merge(mel_nlp, variants[, c("chr", "pos", "codon_start_pos", "classification")], by=c("chr", "pos"), all.x=T)

nrow(variants)
variants <- unique(variants, by = c("chr", "pos"))
nrow(variants)

# nrow(sim_nlp[!is.na(classification)])
# nrow(mel_nlp[!is.na(classification)])

nrow(mel_nlp) 
nrow(sim_nlp) 



mel_dt <- merge(mel_nlp[, .(chr, pos, variant, nLocales_poly, global_af, poly_af, poly_maf)], variants, by = c("chr", "pos"), all.x=TRUE)
sim_dt <- merge(sim_nlp[, .(chr, pos, variant, nLocales_poly, global_af, poly_af, poly_maf)], variants, by = c("chr", "pos"), all.x=TRUE)

af_threshold <- 0.10
maf_label <- af_threshold * 100

mel_dt <- mel_dt[poly_af > af_threshold & poly_af < (1 - af_threshold)]
sim_dt <- sim_dt[poly_af > af_threshold & poly_af < (1 - af_threshold)]

mel_dt[, classification := NULL]
sim_dt[, classification := NULL]


# mel_only_cols <- c(
#   "ref_mel", "alt_mel",
#   "nt_change_mel",
#   "codon_ref_mel", "codon_alt_mel",
#   "aa_ref_mel", "aa_alt_mel",
#   "af_mel",
#   "effect_mel",
#   "gene_mel",
#   "gene_id_mel",
#   "aa_pos_mel"
# )

# sim_only_cols <- c(
#   "ref_sim", "alt_sim",
#   "nt_change_sim",
#   "codon_ref_sim", "codon_alt_sim",
#   "aa_ref_sim", "aa_alt_sim",
#   "af_sim",
#   "effect_sim",
#   "gene_sim",
#   "gene_id_sim",
#   "aa_pos_sim"
# )



suffix_cols <- c("variant", "nLocales_poly", "global_af", "poly_af", "poly_maf")

mel_dt[, (suffix_cols) := lapply(.SD, identity), .SDcols = suffix_cols]
sim_dt[, (suffix_cols) := lapply(.SD, identity), .SDcols = suffix_cols]
setnames(mel_dt, suffix_cols, paste0(suffix_cols, "_mel"))
setnames(sim_dt, suffix_cols, paste0(suffix_cols, "_sim"))

shared_cols <- intersect(
  setdiff(names(mel_dt), c("chr", "pos", paste0(suffix_cols, "_mel"))),
  setdiff(names(sim_dt), c("chr", "pos", paste0(suffix_cols, "_sim")))
)

mel_dt <- mel_dt[, c("chr", "pos", paste0(suffix_cols, "_mel"), shared_cols), with=FALSE]
sim_dt <- sim_dt[, c("chr", "pos", paste0(suffix_cols, "_sim"), shared_cols), with=FALSE]


# id_cols <- setdiff(names(mel_dt), c("chr","pos"))
# shared_cols <- intersect(
#     setdiff(names(mel_dt), c("chr","pos", paste0(suffix_cols, "_mel"))),
#     setdiff(names(sim_dt), c("chr","pos", paste0(suffix_cols, "_sim")))
# )

# mel_dt <- mel_dt[, c("chr","pos", paste0(suffix_cols, "_mel"), shared_cols), with=FALSE]
# sim_dt <- sim_dt[, c("chr","pos", paste0(suffix_cols, "_sim"), shared_cols), with=FALSE]

mel_dt <- mel_dt[!is.na(ref_mel)]
sim_dt <- sim_dt[!is.na(ref_sim)]

sim_only_cols <- names(mel_dt)[endsWith(names(mel_dt), "_sim")]
mel_only_cols <- names(sim_dt)[endsWith(names(sim_dt), "_mel")]

mel_dt <- mel_dt[, (sim_only_cols) := NULL]
sim_dt <- sim_dt[, (mel_only_cols) := NULL]

# mel_dt[, PostMode := NULL]
# mel_dt[, PostMedian := NULL]

# sim_dt[, PostMode := NULL]
# sim_dt[, PostMedian := NULL]

identical(
  mel_dt[, .(chr,pos,codon_start_pos)],
  sim_dt[, .(chr,pos,codon_start_pos)]
)

voi <- merge(
    mel_dt, 
    sim_dt,
    by = c("chr", "pos"),
    all =TRUE
)

voi[
  !is.na(codon_start_pos.x) &
  !is.na(codon_start_pos.y) &
  codon_start_pos.x != codon_start_pos.y,
  .N
]

voi[, codon_start_pos := codon_start_pos.x]
voi[, codon_start_pos.x := NULL]
voi[, codon_start_pos.y := NULL]

names(voi)

file_name <- paste0("/project/berglandlab/anjali/drosophila_polymorphism/classification/MAF", 
  maf_label,
  "filter/subset_qualVar_ofInterest_MAF", 
  maf_label, ".rds")
saveRDS(voi, file_name)


 
# old voi, made from master candidate file, applied MAF 5% to it
# maf5_voi <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/old_tables/subset_qualVar_ofInterest_MAF5_06-18-2026.rds")

# next - run classification_final
## no need to run filteredRDS_varOfInterest.R because ran before maf filter (input rds should be output of filteredRDS_varOfInterest)

# saveRDS(voi, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_MAF5.rds")
# saveRDS(voi_candidates, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/currentFiles/subset_qualVar_ofInterest_classed_geva_MAF5_ABFGOPXY.rds")


#########################################################################################


# background table, not classed, no candidate polymorphisms selected yet
variants <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/merged_tables/quality/all_quality_variants_merge_unfilt.rds")
nrow(variants)

variants <- unique(variants, by = c("chr", "pos"))
nrow(variants)

nrow(mel_nlp) 
nrow(sim_nlp) 

mel_dt <- merge(mel_nlp[, .(chr, pos, variant, nLocales_poly, global_af, poly_af, poly_maf)], variants, by = c("chr", "pos"), all.x=TRUE)
sim_dt <- merge(sim_nlp[, .(chr, pos, variant, nLocales_poly, global_af, poly_af, poly_maf)], variants, by = c("chr", "pos"), all.x=TRUE)

mel_dt <- mel_dt[poly_af > 0.05 & poly_af < 0.95]
sim_dt <- sim_dt[poly_af > 0.05 & poly_af < 0.95]


suffix_cols <- c("variant", "nLocales_poly", "global_af", "poly_af", "poly_maf")

mel_dt[, (suffix_cols) := lapply(.SD, identity), .SDcols = suffix_cols]
sim_dt[, (suffix_cols) := lapply(.SD, identity), .SDcols = suffix_cols]
setnames(mel_dt, suffix_cols, paste0(suffix_cols, "_mel"))
setnames(sim_dt, suffix_cols, paste0(suffix_cols, "_sim"))

shared_cols <- intersect(
  setdiff(names(mel_dt), c("chr", "pos", paste0(suffix_cols, "_mel"))),
  setdiff(names(sim_dt), c("chr", "pos", paste0(suffix_cols, "_sim")))
)

mel_dt <- mel_dt[, c("chr", "pos", paste0(suffix_cols, "_mel"), shared_cols), with=FALSE]
sim_dt <- sim_dt[, c("chr", "pos", paste0(suffix_cols, "_sim"), shared_cols), with=FALSE]


mel_dt <- mel_dt[!is.na(ref_mel)]
sim_dt <- sim_dt[!is.na(ref_sim)]

sim_only_cols <- names(mel_dt)[endsWith(names(mel_dt), "_sim")]
mel_only_cols <- names(sim_dt)[endsWith(names(sim_dt), "_mel")]


mel_dt <- mel_dt[, (sim_only_cols) := NULL]
sim_dt <- sim_dt[, (mel_only_cols) := NULL]


background <- merge(
    mel_dt, 
    sim_dt,
    by = c("chr", "pos"),
    all =TRUE
)

names(background)
saveRDS(voi, "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/merged_tables/quality/all_quality_variants_MAF5_merge_unfilt.rds")

#########################################################################################
