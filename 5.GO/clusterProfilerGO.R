
library(clusterProfiler)
library(org.Dm.eg.db)
library(enrichplot)
library(DOSE)
library(ggplot2)

### SET UP (similar to gowinda, but gene-level not snp-level)

### building background from snps with 1 polymorphic site per codon 
# includes single-species and both-species polymorphic codons
rds <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/noMAFfilter/all_quality_variants_clean.rds")

rds[, `:=`(
  codon_ref_use = fifelse(!is.na(codon_ref_mel), codon_ref_mel, codon_ref_sim),
  nt_ref_use = fifelse(!is.na(ref_mel), ref_mel, ref_sim)
)]

rds[, strand := fifelse(
  nt_ref_use == toupper(substr(codon_ref_use, regexpr("[A-Z]", codon_ref_use), 1)), "forward", "reverse"
)]

rds[strand == "forward", codon_start_pos := pos - (regexpr("[A-Z]", codon_ref_use) - 1)]
rds[strand == "reverse", codon_start_pos := pos + (regexpr("[A-Z]", codon_ref_use) - 1)]
rds[, c("codon_ref_use", "nt_ref_use", "strand") := NULL]
# setkey(filtered_dt, chr, pos)

mel_only <- rds[!is.na(ref_mel)]
sim_only <- rds[!is.na(ref_sim)]

mel_codon_count <- mel_only[, .(n_mel = .N), by = .(chr, codon_start_pos)]
sim_codon_count <- sim_only[, .(n_sim = .N), by = .(chr, codon_start_pos)]

bg_mel <- merge(mel_only, mel_codon_count[n_mel==1], by = c("chr", "codon_start_pos"))[, n_mel := NULL]
bg_sim <- merge(sim_only, sim_codon_count[n_sim==1], by = c("chr", "codon_start_pos"))[, n_sim := NULL]

# includes species-specific codons for the question: 
    # are TSP/conv enriched among all polymorphisms (where poly defined as 1 poly site per codon)
bg_SpeciesSpecific <- rbindlist(list(bg_mel, bg_sim), use.names=T)

# filter for codons where each species is polymorphic at a single site per codon
    # are TSP/conv enriched among shared polymorphisms for specific GO categories
bg_SharedOnly <- merge(mel_codon_count[n_mel == 1], sim_codon_count[n_sim == 1],
  by = c("chr", "codon_start_pos")
)
bg_SharedOnly <- merge(rds, bg_SharedOnly[, .(chr, codon_start_pos)], by = c("chr", "codon_start_pos"))

bg_species_genes <- unique(na.omit(bg_SpeciesSpecific$gene_id_mel))
bg_shared_genes  <- unique(na.omit(bg_SharedOnly$gene_id_mel))

saveRDS(bg_species_genes, "/scratch/ejy4bu/drosophila/GO/clusterProfiler/bg_speciesSpecific_genes.rds")

saveRDS(bg_shared_genes, "/scratch/ejy4bu/drosophila/GO/clusterProfiler/bg_sharedOnly_genes.rds")

#### make candidate files, varying by MAF thresholds 

### Set up: load nlp files
# sim file
load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_simulans.17_06_2026.nlpTable.Rdata")
# load("/project/berglandlab/multispecies_endemism/data/collectiveAnalysis_version3/Drosophila_simulans.10_03_2026.nlpTable.paramask.genmap.busco.Rdata")
### this sim file is in simulans genospace, wrong coordinates
sim_nlp <- nlp 
rm(nlp)

# mel file
# load("/project/berglandlab/anjali/drosophila_polymorphism/data_files/nlp/Drosophila_melanogaster.17_06_2026.nlpTable.Rdata")
load("/project/berglandlab/multispecies_endemism/data/collectiveAnalysis_version3/Drosophila_melanogaster.10_03_2026.nlpTable.paramask.genmap.busco.Rdata")
mel_nlp <- nlp
rm(nlp)   

masterCandidates <- readRDS("/scratch/ejy4bu/drosophila/gowinda/candidateFiles/masterCandidateFile.rds")

# all variant table, no maf filter yet:
sim_nlp <- merge(sim_nlp, masterCandidates[, c("chr", "pos", "codon_start_pos", "classification")], by=c("chr", "pos"), all.x=T)
mel_nlp <- merge(mel_nlp, masterCandidates[, c("chr", "pos", "codon_start_pos", "classification")], by=c("chr", "pos"), all.x=T)

nrow(masterCandidates)

mel_dt <- merge(mel_nlp[, .(chr, pos, variant, nLocales_poly, global_af, poly_af, poly_maf)], masterCandidates, by = c("chr", "pos"), all.x=FALSE)
sim_dt <- merge(sim_nlp[, .(chr, pos, variant, nLocales_poly, global_af, poly_af, poly_maf)], masterCandidates, by = c("chr", "pos"), all.x=FALSE)


getCandidates_filterByMAF <- function(maf, setting) {
    ### testing: params:
    # maf = 0.02
    # setting = "poly"

    if ((maf <= 0) || (setting%in%c("global", "poly") == FALSE)) stop("invalid params")

    classesOfInterest <- c("A", "B", "F", "G", "O", "P", "X", "Y")
    mel_dt <- copy(mel_dt)
    sim_dt <- copy(sim_dt)

    maf_label <- maf * 100

    if (setting=="global") {
        mel_dt <- mel_dt[global_af > maf & global_af < (1 - maf)]
        sim_dt <- sim_dt[global_af > maf & global_af < (1 - maf)]
    } else {
        mel_dt <- mel_dt[poly_af > maf & poly_af < (1 - maf)]
        sim_dt <- sim_dt[poly_af > maf & poly_af < (1 - maf)]
    }

    mel_dt <- mel_dt[!is.na(ref_mel)]
    sim_dt <- sim_dt[!is.na(ref_sim)]

    mel_dt[, mel_present := T]
    sim_dt[, sim_present := T]


    shared_dt <- merge (
        mel_dt[, .(chr, pos, codon_start_pos, classification)], 
        sim_dt[, .(chr, pos, codon_start_pos, classification)], 
        by = c("chr", "codon_start_pos", "classification"), 
        suffixes = c("_mel", "_sim"),
        all = FALSE
    )

    mel_pos <- shared_dt[, .(chr, pos = pos_mel, classification)]
    sim_pos <- shared_dt[, .(chr, pos = pos_sim,classification)]

    all_pos <- rbind(mel_pos, sim_pos)
    unique_chr_pos <- unique(all_pos[!is.na(chr) & !is.na(pos)])
    setorder(unique_chr_pos, chr, pos)

    ### save candidates
    dir = paste0("/scratch/ejy4bu/drosophila/gowinda/candidateFiles/MAF", maf_label, "filter_", setting, "AF/") 

    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
    }
    fwrite(unique_chr_pos[classification%in%classesOfInterest, .(chr, pos)], 
        paste0(dir, "candidate_chrpos_ABFGOPXY_", maf_label, "_", setting, "AF.txt"), sep="\t", col.names=FALSE)
    fwrite(unique_chr_pos[classification%in%c("A", "B"), .(chr, pos)], 
        paste0(dir, "candidate_chrpos_AB_", maf_label, "_", setting, "AF.txt"), sep="\t", col.names=FALSE)
    fwrite(unique_chr_pos[classification%in%c("F", "G", "O", "P", "X", "Y"), .(chr, pos)], 
        paste0(dir, "candidate_chrpos_FGOPXY_", maf_label, "_", setting, "AF.txt"), sep="\t", col.names=FALSE)
    fwrite(unique_chr_pos[classification%in%c("X", "Y"), .(chr, pos)], 
        paste0(dir, "candidate_chrpos_XY_", maf_label, "_", setting, "AF.txt"), sep="\t", col.names=FALSE)
}

maf_inputs <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.49)

for (m in maf_inputs) {
    getCandidates_filterByMAF(maf = m, "poly")
    getCandidates_filterByMAF(maf = m, "global")
}
