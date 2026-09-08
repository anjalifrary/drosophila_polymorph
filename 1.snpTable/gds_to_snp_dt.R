library(SeqArray)
library(data.table)
library(foreach)
library(doMC)

######################################################################
# ### Pool seq files ###
# out_dir <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/"
# out_rds <- paste0(out_dir, "sim_eff_snp_dt.rds")
# if(!file.exists(out_rds)) file.create(out_rds)
# out_csv <- paste0(out_dir, "sim_eff_snp_dt.csv")
# if(!file.exists(out_csv)) file.create(out_csv)

# # load gds file
# # mel_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.eff.gds"
# # gds_file <- seqOpen(mel_file)
# sim_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.eff.gds"
# gds_file <- seqOpen(sim_file)

#######################################################################

### INBRED ###

# ### sim - signor -  inbred 
# out_dir <- "/scratch/ejy4bu/drosophila/inbred/snpDT/"
# out_rds <- paste0(out_dir, "dsim3.signor.snp_dt.rds")
# gds_file <- seqOpen("/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.goodSamps.goodSites.gds")

### mel - DGRP2 - inbred 
out_dir <- "/scratch/ejy4bu/drosophila/inbred/snpDT/"
filtered_rds <- paste0(out_dir, "DGRP2.source_BCM-HGSC.dm6.snp_dt_SynMissense.rds")
full_rds <- paste0(out_dir, "DGRP2.source_BCM-HGSC.dm6.snp_dt_allEffects.rds")
gds_file <- seqOpen("/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.goodSamps.goodSites.gds")



######################################################################

### clean script 9/8/26 ###

filter_effects <- c("synonymous_variant", "missense_variant")

seqResetFilter(gds_file)

dt <- data.table(
    chr = seqGetData(gds_file, "chromosome"), 
    pos = seqGetData(gds_file, "position"),
    nAlleles = seqGetData(gds_file, "$num_allele"),
    id = seqGetData(gds_file, "variant.id"),
    af = seqGetData(gds_file, "annotation/info/AF"),
    maf = seqGetData(gds_file, "annotation/info/MAF"),
    n_samps = seqGetData(gds_file, "annotation/info/NS")
)
dt[, count_records := .N, by = .(chr, pos)]

nrow(dt) # signor: 2938460
summary(dt$af)
summary(dt$maf)

nrow(dt[maf>0,])
nrow(dt[af>0 & af<1, ])
# there are 2630 sites in signor dt that are fixed for ALT allele 

biallelic_dt <- dt[count_records==1 & nAlleles == 2, ] # gets 1 record per (chr, pos) where records have 2 alleles each
nrow(biallelic_dt) # signor: 2830779

biallelic_dt <- biallelic_dt[maf>0, ] # removed fake biallelic records 
nrow(biallelic_dt) # signor: 2828149

seqSetFilter(gds_file, variant.id = biallelic_dt$id)

variant_ids <- biallelic_dt$id

bin_size <- length(variant_ids) # test on 100 variants first 
bins <- split(seq_along(variant_ids), ceiling(seq_along(variant_ids) / bin_size))
n_bins <- length(bins)

# result <- foreach(i = seq_along(bins), .combine = rbind) %do% {   # for binning tables... 
    i = 1

    message("Bin ", i , "/", n_bins)
    idx <- bins[[i]]
    bin_ids <- variant_ids[idx]
    seqSetFilter(gds_file, variant.id = bin_ids)

    # message("getting alleles")
    alleles_all <- seqGetData(gds_file, "allele")
    allele_split <- tstrsplit(alleles_all, ",")

    snp.dt1 <- data.table(
        variant.id = bin_ids,
        chr        = biallelic_dt$chr[idx],
        pos        = biallelic_dt$pos[idx],
        ref        = allele_split[[1]],
        alt        = allele_split[[2]],
        af         = biallelic_dt$af[idx],
        maf        = biallelic_dt$maf[idx],
        n_samps    = biallelic_dt$n_samps[idx]
    )

    # EFF, ., String, Predicted effects for this variant.
    # Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| 
    # (5) Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | 
    # (10) Exon_Rank  | Genotype [ | ERRORS | WARNINGS ] )'
    
    eff_all <- seqGetData(gds_file, "annotation/info/EFF")
    annotated_ids <- seqGetData(gds_file, "variant.id") # keep only annotated variants
    
    eff_dt <- data.table(
        variant.id = rep(annotated_ids, times = eff_all$length),
        eff = eff_all$data
    )
    # keep highest priority snpEff annotation:
    eff_dt[, effect_order := seq_len(.N), by = variant.id]
    eff_row1 <- eff_dt[effect_order == 1, ]

    # extract data
    eff_row1[, effect := sub("\\(.*$", "", eff)]

    eff_row1 <- eff_row1[effect%in%(filter_effects)]

    eff_row1[, eff_contents := sub("^[^(]*\\((.*)\\)$", "\\1", eff)]
    eff_split <- tstrsplit(eff_row1$eff_contents, "\\|")

    eff_row1[, impact           := eff_split[[1]]]  # high/moderate/low/modifier
    eff_row1[, functional_class := eff_split[[2]]]  
    eff_row1[, codon_change     := eff_split[[3]]]  # codon change. format : aAt/aCt
    eff_row1[, aa_change        := eff_split[[4]]]  # aa change. p.Gln5His
    eff_row1[, aa_length        := eff_split[[5]]]  # aa length as integer
    eff_row1[, gene             := eff_split[[6]]]  # gene, format CG11023
    eff_row1[, biotype          := eff_split[[7]]]  # biotype like protein coding (should all be protein coding)
    eff_row1[, gene_coding      := eff_split[[8]]]  # gene coding "CODING"
    eff_row1[, transcript_id    := eff_split[[9]]]  # transcript like FBtr000000000
    eff_row1[, exon_rank        := eff_split[[10]]] # integer describing position of exon in transcript
    eff_row1[, genotype         := eff_split[[11]]] # nucleotide

    eff_row1[, eff_contents := NULL]
    eff_row1[, eff := NULL]
    eff_row1[, effect_order := NULL]

    snp_dt <- merge(snp.dt1, eff_row1, by = "variant.id", all.x=T)
    filtered_dt <- snp_dt[effect%in%(filter_effects)]

    saveRDS(snp_dt, full_rds)
    saveRDS(filtered_dt, filtered_rds)
    seqClose(gds_file)


# }

####################################################################

### old script - obsolete 

# filter_effects <- c("synonymous_variant", "missense_variant")

# ### loading both dt tables
# # Filter for biallelic variants
# build_snp_dt <- function(gds) {
#     seqResetFilter(gds)
#     dt <- data.table(
#         chr      = seqGetData(gds, "chromosome"),
#         pos      = seqGetData(gds, "position"),
#         nAlleles = seqGetData(gds, "$num_allele"),
#         id       = seqGetData(gds, "variant.id")
#     )
#     dt <- dt[nAlleles == 2]
#     seqSetFilter(gds, variant.id = dt$id)
#     dt[, af := seqGetData(gds, "annotation/info/AF")]
#     message(nrow(dt), " biallelic variants")
#     dt
# }


# ### extract annotations for each variant
# build_species_dt <- function(gds, snp_dt, bin_size=10000){
#     snp_id <- snp_dt$id
#     bins <- split(seq_along(snp_id), ceiling(seq_along(snp_id) / bin_size))
#     n_bins <- length(bins)

#     result <- foreach(i = seq_along(bins), .combine = rbind) %do% {
#         message("Bin ", i , "/", n_bins)
#         idx <- bins[[i]]
#         bin_id <- snp_id[idx]
#         seqSetFilter(gds, variant.id = bin_id)

#         # message("getting alleles")
#         alleles_all <- seqGetData(gds, "allele")
#         allele_split <- tstrsplit(alleles_all, ",")

#         snp.dt1 <- data.table(
#             variant.id = bin_id,
#             chr        = snp_dt$chr[idx],
#             pos        = snp_dt$pos[idx],
#             ref        = allele_split[[1]],
#             alt        = allele_split[[2]],
#             af         = snp_dt$af[idx])

#         # message("getting annotations")
#         ann_all <- seqGetData(gds, "annotation/info/ANN")
#         # message("ann_all length field: ", length(ann_all$length))
#         # message("bin_id length: ", length(bin_id))

#         annotated_ids <- seqGetData(gds, "variant.id")  # filter out if no annotation

#         ann_dt <- data.table( variant.id = rep(annotated_ids, times=ann_all$length), ann = ann_all$data)
#         # add effect order 
#         ann_dt[, effect_order := seq_len(.N), by = variant.id]

#         ann_split <- tstrsplit(ann_dt$ann, "\\|")

#         ann_dt[,effect := ann_split[[2]]]           #class of annotation (e.g. upstream_gene_variant)
#         ann_dt[, impact := ann_split[[3]]]          # high/moderate/low/modifier
#         ann_dt[, gene := ann_split[[4]]]            # gene name
#         ann_dt[, gene_id := ann_split[[5]]]         # flybase gene id
#         ann_dt[, feature_type := ann_split[[6]]]    # e.g. transcript
#         ann_dt[, transcript_id := ann_split[[7]]]   #
#         ann_dt[, biotype := ann_split[[8]]]         # e.g. protein-coding
#         ann_dt[, in_exon_pos := ann_split[[9]]]     # intron or exon position
#         ann_dt[, nt_change := ann_split[[10]]]      # nucleotide change & position (c.-1427T>A)
#         ann_dt[, nt_pos := ann_split[[13]]]         # amino acid position within the protein

#         ann_dt[, aa_change := ann_split[[11]]]      # amino acid change

#         ann_dt[, ann := NULL]  # drop the raw string, keep parsed columns

#         eff_all <- seqGetData(gds, "annotation/info/EFF")
#         eff_dt <- data.table( variant.id = rep(annotated_ids, times=eff_all$length), eff = eff_all$data)
#         eff_split <- tstrsplit(eff_dt$eff, "\\|")
#         eff_dt[, codon_change := eff_split[[3]]]    # codon change. format : aAt/aCt

#         # keep first eff per variant
#         eff_dt_first <- eff_dt[, .SD[1], by = variant.id]

#         # merge EFF with ANN
#         ann_dt <- merge(ann_dt, eff_dt_first[, .(variant.id, codon_change)], by="variant.id", all.x=TRUE)

#         aa_consistent_bin <- ann_dt[aa_change != "", .(
#             num_transcripts = .N,
#             num_unique_aa = uniqueN(aa_change),
#             consistent = uniqueN(aa_change) == 1
#         ), by = variant.id]
        
#         ## merge binned table
#         bin_table <- merge(snp.dt1, ann_dt[, .(variant.id, effect_order, effect, impact, gene, gene_id, 
#         feature_type, transcript_id, biotype, in_exon_pos, nt_change, nt_pos, aa_change, codon_change)], by = "variant.id")

#         ## keep only first variant (by effect_order column)
#         bin_table <- bin_table[effect_order==1]
#         # ann_canonical <- ann_dt[, effect_order==1]

#         rm(ann_all, ann_dt, ann_split, snp.dt1)
#         gc()
#         bin_table
#     }
#     return(result)
# }

# snp_dt <- build_snp_dt(gds_file)
# species_table <- build_species_dt(gds_file, snp_dt)

# # filter for synonymous or missense 
# species_table <- species_table[effect %in% filter_effects] 
# message("filtered variants: kept ", filter_effects)

# # check that all amino acid polymorphisms are the same even with different transcripts:
# aa_consistent <- species_table[aa_change != "", .(
#     num_transcripts = .N,
#     num_unique_aa = uniqueN(aa_change),
#     consistent = uniqueN(aa_change)==1
# ), by = variant.id]
# message("variants with consistent aa_change: ", sum(aa_consistent$consistent))
# message("variants with inconsistent aa_change: ", sum(!aa_consistent$consistent))
# aa_consistent[consistent == FALSE][1:10] # view first 10 inconsistent variants

# saveRDS(species_table, out_rds)
# message("variants: ", nrow(species_table), "\nsaved to: ", out_rds)

# # subset_table <- species_table[1:500, ]
# # fwrite(subset_table, out_csv)
# # message("saved first 500 rows to csv at ", out_csv)
