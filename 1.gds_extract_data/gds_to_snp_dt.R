library(SeqArray)
library(data.table)
library(foreach)
library(doMC)

out_dir <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_datatables/test_files/"
out_csv <- paste0(out_dir, "mel_snp_dt.csv")
out_rds <- paste0(out_dir, "mel_snp_dt.rds")
if(!file.exists(out_csv)) file.create(out_csv)
if(!file.exists(out_rds)) file.create(out_rds)


# load melanogaster gds file
mel_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.gds"
# sim_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds"
gds_file <- seqOpen(mel_file)
# gds_file <- seqOpen(sim_file)

# genofile <- sim_file

# genofile <- seqOpen(genofile)
# genofile

# # load simulans gds file
# sim_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds"
# sim_gds <- seqOpen(sim_file)
# sim_gds 

filter_effects <- c("synonymous_variant", "missense_variant")

### loading both dt tables
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
build_species_dt <- function(gds, snp_dt, bin_size=2000){
    snp_id <- snp_dt$id
    bins <- split(seq_along(snp_id), ceiling(seq_along(snp_id) / bin_size))
    n_bins <- length(bins)

    result <- foreach(i = seq_along(bins), .combine = rbind) %do% {
        message("Bin ", i , "/", n_bins)
        idx <- bins[[i]]
        bin_id <- snp_id[idx]
        seqSetFilter(gds, variant.id = bin_id)

        message("getting alleles")
        alleles_all <- seqGetData(gds, "allele")
        allele_split <- tstrsplit(alleles_all, ",")

        snp.dt1 <- data.table(
            variant.id = bin_id,
            chr        = snp_dt$chr[idx],
            pos        = snp_dt$pos[idx],
            ref        = allele_split[[1]],
            alt        = allele_split[[2]],
            af         = snp_dt$af[idx])

        message("getting annotations")
        ann_all <- seqGetData(gds, "annotation/info/ANN")
        message("ann_all length field: ", length(ann_all$length))
        message("bin_id length: ", length(bin_id))

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
        ann_dt[, biotype := ann_split[[8]]]         #e.g. protein-coding
        ann_dt[, in_exon := ann_split[[9]]]         # intron or exon
        ann_dt[, nt_change := ann_split[[10]]]      # nucleotide change & position (c.-1427T>A)
        ann_dt[, aa_change := ann_split[[11]]]      # amino acid change
        ann_dt[, aa_pos := ann_split[[13]]]         # amino acid position within the protein

        ann_dt[, ann := NULL]  # drop the raw string, keep parsed columns

        ## merge binned table
        bin_table <- merge(snp.dt1, ann_dt[, .(variant.id, effect, impact, gene, gene_id, 
        feature_type, transcript_id, biotype, in_exon, nt_change, aa_change, aa_pos)], by = "variant.id")

        ## keep only first variant (by effect_order column)


        rm(ann_all, ann_dt, ann_split, snp.dt1)
        gc()
        bin_table
    }
    return(result)
}

#### test
mel_snp_dt <- build_snp_dt(gds_file)
mel_snp_dt_test <- mel_snp_dt[1:1000]
mel_table <- build_species_dt(gds_file, mel_snp_dt_test)

# check that all amino acid polymorphisms are the same even with different transcripts:
aa_consistent <- mel_table[aa_change != "", .(
    num_transcripts = .N,
    num_unique_aa = uniqueN(aa_change),
    consistent = uniqueN(aa_change)==1
), by = variant.id]
message("variants with consistent aa_change: ", sum(aa_consistent$consistent))
message("variants with inconsistent aa_change: ", sum(!aa_consistent$consistent))
aa_consistent[consistent == FALSE][1:10] # view first 10 inconsistent variants

# filter for synonymous or missense 
mel_table <- mel_table[effect %in% filter_effects] 

# mel_table <- mel_table[order(variant.id, biotype = "protein_coding"), .SD[1], by = variant.id] # compress to first variant 
message("variants: ", nrow(mel_table), "\nsaved to: ", out_csv)
saveRDS(mel_table, out_rds)
fwrite(mel_table, out_csv)

# #### build table  ( commented out for tested subset )
# message("----building table-----")
# species_dt <- build_snp_dt(gds_file)
# species_table <- build_species_dt(gds_file, species_dt)
# message("coding variants: ", nrow(species_table))
# saveRDS(species_table, paste0(out_rds))

#### shared table code --- OLD --- ########################################################

# # warnings()
# # mel_snp_dt <- build_snp_dt(mel_gds)
# # sim_snp_dt <- build_snp_dt(sim_gds)

# # build shared datatable
# shared <- merge(mel_snp_dt, sim_snp_dt, by = c("chr", "pos"), suffixes = c("_mel", "_sim"), all=T)
# message(nrow(shared), " total variants")

# # shared_test <- shared[1:100]
# shared_table <- get_gds_data(mel_gds, shared, "mel")
# # message("saving csv to ", out_csv)
# # fwrite(shared_table, out_csv)
# message("saving rds to ", out_rds)
# saveRDS(shared_table, out_rds)

# message("complete. ", nrow(shared_table), " variants written.")

#### end of shared table code --- OLD --- #######################################################

### ---------- uncompressed code ----------- ####################################################

# snp.dt <- data.table(
#     chr=seqGetData(genofile, "chromosome"),
#     pos=seqGetData(genofile, "position"),
#     nAlleles=seqGetData(genofile, "$num_allele"),
#     id=seqGetData(genofile,"variant.id"))

# snp.dt <- snp.dt[nAlleles==2] ##subset to sites with 2 alleles
# seqResetFilter(genofile)
# seqSetFilter(genofile, variant.id = snp.dt$id)  # biallelic only

# ### SHARED SITES
# message("finding shared snps")
# shared <- merge(mel_snp_dt, sim_snp_dt, by=c("chr", "pos"), suffixes=c("mel", "sim"))
# message(nrow(shared), " shared positions")

# message("gettting mel alleles")
# seqSetFilter(mel_gds, variant.id = shared$id_mel)
# alleles_all <- seqGetData(mel_gds, "allele")
# allele_split <- tstrsplit(alleles_all, ",")


# snp.dt1 <- data.table(
#     variant.id = snp.dt$id,
#     chr        = snp.dt$chr,
#     pos        = snp.dt$pos,
#     ref        = allele_split[[1]],
#     alt        = allele_split[[2]])

# message("getting annotations")
# ann_all <- seqGetData(genofile, "annotation/info/ANN")

# ann_dt <- data.table( variant.id = rep(snp.dt$id, times=ann_all$length), 
#     ann = ann_all$data)
# ann_split <- tstrsplit(ann_dt$ann, "\\|")

#  $ length: int 15
#  $ data  : chr [1:15] "A|upstream_gene_variant|MODIFIER|CG11023|FBgn0031208|transcript|FBtr0300690|protein_coding||c.-1427T>A|||||1276|" 
# "A|upstream_gene_variant|MODIFIER|CG11023|FBgn0031208|transcript|FBtr0300689|protein_coding||c.-1427T>A|||||1276|" 
# "A|upstream_gene_variant|MODIFIER|CG11023|FBgn0031208|transcript|FBtr0330654|protein_coding||c.-1427T>A|||||1276|" 
# "A|downstream_gene_variant|MODIFIER|l(2)gl|FBgn0002121|transcript|FBtr0078170|protein_coding||c.*4962A>T|||||3586|" ...   
#  - attr(*, "class")= chr "SeqVarDataList"

# ann_dt[,effect := ann_split[[2]]]           #class of annotation (e.g. upstream_gene_variant)
# ann_dt[, impact := ann_split[[3]]]          # high/moderate/low/modifier
# ann_dt[, gene := ann_split[[4]]]            # gene name
# ann_dt[, gene_id := ann_split[[5]]]         # flybase gene id
# ann_dt[, feature_type := ann_split[[6]]]    # e.g. transcript
# ann_dt[, transcript_id := ann_split[[7]]]   #
# ann_dt[, biotype := ann_split[[8]]]         #e.g. protein-coding
# ann_dt[, in_exon := ann_split[[9]]]         # intron or exon
# ann_dt[, nt_change := ann_split[[10]]]      # nucleotide change & position (c.-1427T>A)
# ann_dt[, aa_change := ann_split[[11]]]      # amino acid change
# ann_dt[, aa_pos := ann_split[[13]]]         # amino acid position within the protein

# message("classifying as syn/nonsyn")

# message("intron/exon classification")

# message("collapse to canonical transcript...")
# ann_canonical <- ann_dt[, .SD[1], by = variant.id]

# # join with SNP info
# snp_table <- merge(snp.dt1, ann_canonical[, .(variant.id, effect, 
#     impact, gene, gene_id, feature_type, transcript_id, biotype, in_exon, 
#     nt_change, aa_change, aa_pos)],by = "variant.id")
