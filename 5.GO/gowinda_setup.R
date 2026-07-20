library(data.table)

# requires four files:
    # background snps in .txt file of chr | pos
    # candidate snps - subset
    # GTF annotation file - make sure I'm using the right one
    # GO gene sets

dir <- "/scratch/ejy4bu/drosophila/gowinda/MAF5/new_6-29-26/"

# ### background snps
# rds <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/all_quality_variants_MAF5_clean.rds")
# total_snp <- unique(rds[, .(chr, pos)])
# fwrite(total_snp, paste0(dir, "background_all_snps.txt"), sep="\t", col.names=FALSE)

### DO ONCE

### gtf annotation file
gtf_file <- fread("/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/dmel-all-r6.67.gtf")

### GO file
gaf <- fread("/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/gene_association.fb",
    sep="\t", header=FALSE, fill=TRUE, skip=5)

# Remove comment lines
# gaf <- gaf[!grepl("^!", V1)]
gowinda_go <- gaf[grepl("^FBgn", V2) & grepl("^GO:", V5),
    .(
        TERM  = unique(V5)[1],
        genes = paste(unique(V2), collapse=" ")
    ), by=V5]
    
fwrite(gowinda_go, "/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/flybase_gaf_go.txt",
    sep="\t", col.names=FALSE, quote=FALSE)

# making GO file
library(AnnotationDbi)
library(org.Dm.eg.db)  # Drosophila melanogaster annotation package
library(data.table)

go_to_fb <- AnnotationDbi::select(org.Dm.eg.db,
    keys    = keys(org.Dm.eg.db, keytype = "GO"),
    columns = c("GO", "FLYBASE"),
    keytype = "GO"
)

go_dt <- as.data.table(go_to_fb)[!is.na(FLYBASE)]

gowinda_go <- go_dt[, .(
    TERM  = unique(GO)[1],   # reuse GO ID as descriptor ** FIX ** 
    genes = paste(unique(FLYBASE), collapse = "\t")
), by = GO]

fwrite(gowinda_go, "/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/flybase_go.txt",
    sep = "\t", col.names = FALSE, quote = FALSE)





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

fwrite(bg_SpeciesSpecific[, .(chr, pos)], "/scratch/ejy4bu/drosophila/gowinda/backgroundFiles/bg_speciesSpecific_noMAF.txt", sep="\t", col.names=F)
fwrite(bg_SharedOnly[, .(chr, pos)], "/scratch/ejy4bu/drosophila/gowinda/backgroundFiles/bg_sharedOnly_noMAF.txt", sep="\t", col.names=F)

### MASTER CANDIDATES FILE - do once
# all shared polymorphisms (classed)
shared_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/noMAFfilter/subset_qualVar_ofInterest_7-20-2026.rds")
classesOfInterest <- c("A", "B", "F", "G", "O", "P", "X", "Y")
masterCandidates <- shared_dt[classification%in%classesOfInterest, ]
saveRDS(masterCandidates, "/scratch/ejy4bu/drosophila/gowinda/candidateFiles/masterCandidateFile.rds")

### FILTER CANDIDATES FOR MAF FILTER - save gowinda inputs

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


# ### filtered background snps
# voi <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/voi_fromBG_qualVar_ofInterest_MAF5_classed_06-29-2026.rds")
# # candidate_snp <- unique(voi[, .(chr, pos)])
# # fwrite(candidate_snp, paste0(dir, "background_classed_snps.txt"), sep="\t", col.names=FALSE)

### candidate snps by classes
# fwrite(unique(voi[classification=="A", .(chr, pos)]), paste0(dir, "candidate_snp_A.txt"), sep="\t", col.names=FALSE)
# fwrite(unique(voi[classification=="B", .(chr, pos)]), paste0(dir, "candidate_snp_B.txt"), sep="\t", col.names=FALSE)
# fwrite(unique(voi[classification%in%c("A", "B"), .(chr, pos)]), paste0(dir, "candidate_snp_AB.txt"), sep="\t", col.names=FALSE)
# fwrite(unique(voi[classification%in%c("F", "G", "O", "P"), .(chr, pos)]), paste0(dir, "candidate_snp_FGOP.txt"), sep="\t", col.names=FALSE)
# fwrite(unique(voi[classification%in%c("X", "Y"), .(chr, pos)]), paste0(dir, "candidate_snp_XY.txt"), sep="\t", col.names=FALSE)
# fwrite(unique(voi[classification%in%c("F", "G", "O", "P", "X", "Y"), .(chr, pos)]), paste0(dir, "candidate_snp_FGOPXY.txt"), sep="\t", col.names=FALSE)
# fwrite(unique(voi[classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y"), .(chr, pos)]), paste0(dir, "candidate_snp_ABFGOPXY.txt"), sep="\t", col.names=FALSE)


# java -Xmx8g -jar Gowinda.jar \
#   --snp-file total_snps.txt \
#   --candidate-snp-file candidate_A.txt \
#   --gene-set-file flybase_go.txt \
#   --annotation-file dmel.gtf \
#   --simulations 1000000 \
#   --gene-definition gene \
#   --threads 16 \
#   --mode gene \
#   --output-file A_gowinda.txt


# for regulatory variation
# --gene-definition updownstream2000

# can set mode to snp instead of gene... should I?