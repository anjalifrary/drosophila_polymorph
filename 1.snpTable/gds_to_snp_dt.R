library(SeqArray)
library(data.table)
library(foreach)
library(doMC)

######################################################################
# ### Pool seq files ###

# mel: 
out_dir <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/mel/"
full_rds <- paste0(out_dir, "dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.snp_dt_allEffects.10_2026.rds")
filtered_rds <- paste0(out_dir, "dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.snp_dt_SynMissense.10_2026.rds")
gds_file <- seqOpen("/scratch/ejy4bu/drosophila/DEST_remake/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.snpGap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.gds")

# sim:
out_dir <- "/scratch/ejy4bu/drosophila/DEST_remake/snpDT/sim/"
full_rds <- paste0(out_dir, "dest.sim.all.SNAPE.001.50.20Nov2025_sim.snp_dt_allEffects.10_2026.rds")
filtered_rds <- paste0(out_dir, "dest.sim.all.SNAPE.001.50.20Nov2025_sim.snp_dt_SynMissense.10_2026.rds")
gds_file <- seqOpen("/scratch/ejy4bu/drosophila/DEST_remake/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.norm.snpGap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.merged.biallelic.dm6.sorted.m2.M2.gds")


#######################################################################

### INBRED ###

# ### sim - signor -  inbred 
# out_dir <- "/scratch/ejy4bu/drosophila/inbred/snpDT/"
# filtered_rds <- paste0(out_dir, "dsim3.signor.snp_dt_SynMissense.rds")
# full_rds <- paste0(out_dir, "dsim3.signor.snp_dt_allEffects.rds")
# gds_file <- seqOpen("/scratch/ejy4bu/drosophila/inbred/gds_files/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.gds")

# ### mel - DGRP2 - inbred 
# out_dir <- "/scratch/ejy4bu/drosophila/inbred/snpDT/"
# filtered_rds <- paste0(out_dir, "DGRP2.source_BCM-HGSC.dm6.snp_dt_SynMissense.rds")
# full_rds <- paste0(out_dir, "DGRP2.source_BCM-HGSC.dm6.snp_dt_allEffects.rds")
# gds_file <- seqOpen("/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.goodSamps.goodSites.gds")


######################################################################

### clean script 9/8/26 ###

filter_effects <- c("synonymous_variant", "missense_variant")

seqResetFilter(gds_file)

# ### this worked on the DGRP (inbred) gds:
# dt <- data.table(
#     chr = seqGetData(gds_file, "chromosome"), 
#     pos = seqGetData(gds_file, "position"),
#     nAlleles = seqGetData(gds_file, "$num_allele"),
#     id = seqGetData(gds_file, "variant.id"),
#     af = seqGetData(gds_file, "annotation/info/AF"),
#     maf = seqGetData(gds_file, "annotation/info/MAF"),
#     n_samps = seqGetData(gds_file, "annotation/info/NS")
# )

### this worked on the dsim3 signor (inbred) gds:
dt <- data.table(
    chr = seqGetData(gds_file, "chromosome"), 
    pos = seqGetData(gds_file, "position"),
    nAlleles = seqGetData(gds_file, "$num_allele"),
    id = seqGetData(gds_file, "variant.id"),
    af = seqGetData(gds_file, "annotation/info/AF")
    # maf = seqGetData(gds_file, "annotation/info/MAF"),
    # n_samps = seqGetData(gds_file, "annotation/info/NS")
)

dt[, count_records := .N, by = .(chr, pos)]

nrow(dt) 
# dgrp: 2830779 #signor: 3905965
# dest sim: 4427142; dest mel: 27302084
# summary(dt$af)
# summary(dt$maf)

# nrow(dt[maf>0,])
# nrow(dt[af>0 & af<1, ])
# there are 2630 sites in dgrp dt that are fixed for ALT allele 

biallelic_dt <- dt[count_records==1 & nAlleles == 2, ] # gets 1 record per (chr, pos) where records have 2 alleles each
nrow(biallelic_dt) # dgrp: 2830779 # signor: 3567947 # dest sim: 4420061 ; dest mel: 27302084

# biallelic_dt <- biallelic_dt[maf>0, ] # removed fake biallelic records 
seqSetFilter(gds_file, variant.id = biallelic_dt$id)

# now get AF for signor:
# biallelic_dt[, af := (seqGetData(gds_file, "annotation/info/AF"))$data]
# biallelic_dt[, maf := pmin(af, 1-af)]

biallelic_dt <- biallelic_dt[maf >0 & af < 1, ] 
nrow(biallelic_dt) # signor: 3566589; dgrp: 2828149
variant_ids <- biallelic_dt$id

# bin_size <- length(variant_ids) # test on 100 variants first 
bin_size <- 1000000
bins <- split(seq_along(variant_ids), ceiling(seq_along(variant_ids) / bin_size))
n_bins <- length(bins)

result <- foreach(i = seq_along(bins), .combine = rbind) %do% {   # for binning tables... 
    # i = 1

    message("Bin ", i , "/", n_bins)
    idx <- bins[[i]]
    bin_ids <- variant_ids[idx]
    seqSetFilter(gds_file, variant.id = bin_ids)

    # message("getting alleles")
    alleles_all <- seqGetData(gds_file, "allele")
    allele_split <- tstrsplit(alleles_all, ",")

    # ### for signor only:
    # genotypes <- seqGetData(gds_file, "genotype")
    # # extract genotypes for all samples and apply (1, 0, 3=missing... and NAs)
    # n_samps <- apply(genotypes, 3, function(g) {
    #     sum(!is.na(g[1, ]) & !is.na(g[2, ]))
    # })

    # src_chr <- seqGetData(gds_file, "annotation/info/SRC_CHROM")
    # src_pos <- seqGetData(gds_file, "annotation/info/SRC_POS")

    # # src_ref_alt is list with length 2 for ALL snps so extracting pair for each snp :
    # src_ref_alt <- seqGetData(gds_file, "annotation/info/SRC_REF_ALT")
    # ref_src <- src_ref_alt$data[seq(1, length(src_ref_alt$data), by = 2)]
    # alt_src <- src_ref_alt$data[seq(2, length(src_ref_alt$data), by = 2)]

    # # flip = if strand was flipped; swapped = if ref and alt were swapped
    # flip <- seqGetData(gds_file, "annotation/info/FLIP")
    # swap <- seqGetData(gds_file, "annotation/info/SWAP")
    # ### end signor specific code ^

    snp.dt1 <- data.table(
        variant.id = bin_ids,
        chr        = biallelic_dt$chr[idx],
        pos        = biallelic_dt$pos[idx],
        ref        = allele_split[[1]],
        alt        = allele_split[[2]],
        af         = biallelic_dt$af[idx]
        # ,
        # maf        = biallelic_dt$maf[idx]
        # ,
        # n_samps    = biallelic_dt$n_samps
        # ,
        # chr_src    = src_chr,
        # pos_src    = src_pos,
        # ref_src    = ref_src,
        # alt_src    = alt_src,
        # flip       = flip,
        # swap       = swap
    )


    # EFF, ., String, Predicted effects for this variant.
    # Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| 
    # (5) Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | 
    # (10) Exon_Rank  | Genotype [ | ERRORS | WARNINGS ] )'
    
    annotated_ids <- seqGetData(gds_file, "variant.id") # keep only annotated variants

    ann_all <- seqGetData(gds_file, "annotation/info/ANN")
    ann_dt <- data.table(
        variant.id = rep(annotated_ids, times = ann_all$length), 
        ann = ann_all$data
    )
    ann_dt[, effect_order := seq_len(.N), by = variant.id]
    ann_row1 <- ann_dt[effect_order==1, ]
    # ann_row1[, ann_contents := sub("^[^(]*\\((.*)\\)$", "\\1", ann)]
    ann_split <- tstrsplit(ann_row1$ann, "\\|")
    ann_row1[, gene_id_fbgn := ann_split[[5]]]

    eff_all <- seqGetData(gds_file, "annotation/info/EFF")

    eff_dt <- data.table(
        variant.id = rep(annotated_ids, times = eff_all$length),
        eff = eff_all$data
    )  
    # keep highest priority snpEff annotation:
    eff_dt[, effect_order := seq_len(.N), by = variant.id]
    eff_row1 <- eff_dt[effect_order == 1, ]

    # extract data
    eff_row1[, effect := sub("\\(.*$", "", eff)]

    # eff_row1 <- eff_row1[effect%in%(filter_effects)]

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
    eff_row1[, genotype         := eff_split[[11]]] # nucleotide snpEff genotype?

    eff_row1[, eff_contents := NULL]
    eff_row1[, eff := NULL]
    eff_row1[, effect_order := NULL]

    snp_dt <- merge(snp.dt1, eff_row1, by = "variant.id", all.x=T)
    snp_dt <- merge(snp_dt, ann_row1[, .(variant.id, gene_id_fbgn)], by = "variant.id", all.x=T)
    # filtered_dt <- snp_dt[effect%in%(filter_effects)]

    bin_table <- merge(snp.dt1, eff_row1, by="variant.id")
    bin_table <- merge(bin_table, ann_row1[, .(variant.id, gene_id_fbgn)], by="variant.id", all.x=T)

    rm(eff_all, ann_all, eff_dt, ann_dt, eff_row1, ann_row1,
        eff_split, ann_split, snp.dt1, alleles_all, allele_split)

    gc()
    bin_table

}
# return(result)
# result <- bin_table

filtered_dt <- result[effect%in%(filter_effects)]

saveRDS(result, full_rds)
saveRDS(filtered_dt, filtered_rds)
seqClose(gds_file)
