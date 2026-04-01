library(SeqArray)
library(data.table)
library(foreach)
library(doMC)

out_dir <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/"
out_rds <- paste0(out_dir, "mel_eff_snp_dt.rds")
if(!file.exists(out_rds)) file.create(out_rds)
out_csv <- paste0(out_dir, "mel_eff_snp_dt.csv")
if(!file.exists(out_csv)) file.create(out_csv)

# load gds file
mel_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.eff.gds"
gds_file <- seqOpen(mel_file)
# sim_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.eff.gds"
# gds_file <- seqOpen(sim_file)

filter_effects <- c("synonymous_variant", "missense_variant")

### loading both dt tables
# Filter for biallelic variants
build_snp_dt <- function(gds) {
    seqResetFilter(gds)
    dt <- data.table(
        chr      = seqGetData(gds, "chromosome"),
        pos      = seqGetData(gds, "position"),
        nAlleles = seqGetData(gds, "$num_allele"),
        id       = seqGetData(gds, "variant.id")
    )
    dt <- dt[nAlleles == 2]
    seqSetFilter(gds, variant.id = dt$id)
    dt[, af := seqGetData(gds, "annotation/info/AF")$data]
    message(nrow(dt), " biallelic variants")
    dt
}


### extract annotations for each variant
build_species_dt <- function(gds, snp_dt, bin_size=10000){
    snp_id <- snp_dt$id
    bins <- split(seq_along(snp_id), ceiling(seq_along(snp_id) / bin_size))
    n_bins <- length(bins)

    result <- foreach(i = seq_along(bins), .combine = rbind) %do% {
        message("Bin ", i , "/", n_bins)
        idx <- bins[[i]]
        bin_id <- snp_id[idx]
        seqSetFilter(gds, variant.id = bin_id)

        # message("getting alleles")
        alleles_all <- seqGetData(gds, "allele")
        allele_split <- tstrsplit(alleles_all, ",")

        snp.dt1 <- data.table(
            variant.id = bin_id,
            chr        = snp_dt$chr[idx],
            pos        = snp_dt$pos[idx],
            ref        = allele_split[[1]],
            alt        = allele_split[[2]],
            af         = snp_dt$af[idx])

        # message("getting annotations")
        ann_all <- seqGetData(gds, "annotation/info/ANN")
        # message("ann_all length field: ", length(ann_all$length))
        # message("bin_id length: ", length(bin_id))

        annotated_ids <- seqGetData(gds, "variant.id")  # filter out if no annotation

        ann_dt <- data.table( variant.id = rep(annotated_ids, times=ann_all$length), ann = ann_all$data)
        # add effect order 
        ann_dt[, effect_order := seq_len(.N), by = variant.id]

        ann_split <- tstrsplit(ann_dt$ann, "\\|")

        ann_dt[,effect := ann_split[[2]]]           #class of annotation (e.g. upstream_gene_variant)
        ann_dt[, impact := ann_split[[3]]]          # high/moderate/low/modifier
        ann_dt[, gene := ann_split[[4]]]            # gene name
        ann_dt[, gene_id := ann_split[[5]]]         # flybase gene id
        ann_dt[, feature_type := ann_split[[6]]]    # e.g. transcript
        ann_dt[, transcript_id := ann_split[[7]]]   #
        ann_dt[, biotype := ann_split[[8]]]         # e.g. protein-coding
        ann_dt[, in_exon_pos := ann_split[[9]]]     # intron or exon position
        ann_dt[, nt_change := ann_split[[10]]]      # nucleotide change & position (c.-1427T>A)
        ann_dt[, nt_pos := ann_split[[13]]]         # amino acid position within the protein

        ann_dt[, aa_change := ann_split[[11]]]      # amino acid change
        # ann_dt[, codon_change := ann_split[[12]]]   # codon that codes for amino acid 

        ann_dt[, ann := NULL]  # drop the raw string, keep parsed columns

        eff_all <- seqGetData(gds, "annotation/info/EFF")
        eff_dt <- data.table( variant.id = rep(annotated_ids, times=eff_all$length), eff = eff_all$data)
        eff_split <- tstrsplit(eff_dt$eff, "\\|")
        eff_dt[, codon_change := eff_split[[3]]]    # codon change. format : aAt/aCt

        # keep first eff per variant
        eff_dt_first <- eff_dt[, .SD[1], by = variant.id]

        # merge EFF with ANN
        ann_dt <- merge(ann_dt, eff_dt_first[, .(variant.id, codon_change)], by="variant.id", all.x=TRUE)

        ## merge binned table
        bin_table <- merge(snp.dt1, ann_dt[, .(variant.id, effect_order, effect, impact, gene, gene_id, 
        feature_type, transcript_id, biotype, in_exon_pos, nt_change, nt_pos, aa_change, codon_change)], by = "variant.id")

        ## keep only first variant (by effect_order column)
        bin_table <- bin_table[effect_order==1]
        # ann_canonical <- ann_dt[, effect_order==1]

        rm(ann_all, ann_dt, ann_split, snp.dt1)
        gc()
        bin_table
    }
    return(result)
}

snp_dt <- build_snp_dt(gds_file)
species_table <- build_species_dt(gds_file, snp_dt)

# filter for synonymous or missense 
species_table <- species_table[effect %in% filter_effects] 
message("filtered variants: kept ", filter_effects)

# # check that all amino acid polymorphisms are the same even with different transcripts:
# aa_consistent <- species_table[aa_sub != "", .(
#     num_transcripts = .N,
#     num_unique_aa = uniqueN(aa_sub),
#     consistent = uniqueN(aa_sub)==1
# ), by = variant.id]
# message("variants with consistent aa_change: ", sum(aa_consistent$consistent))
# message("variants with inconsistent aa_change: ", sum(!aa_consistent$consistent))
# aa_consistent[consistent == FALSE][1:10] # view first 10 inconsistent variants

saveRDS(species_table, out_rds)
message("variants: ", nrow(species_table), "\nsaved to: ", out_rds)

subset_table <- species_table[1:500, ]
fwrite(subset_table, out_csv)
message("saved first 500 rows to csv at ", out_csv)
