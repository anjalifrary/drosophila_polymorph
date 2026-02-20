library(SeqArray)
library(data.table)

out_dir <- "/scratch/ejy4bu/drosophila/gds_analysis/"
out_csv <- paste0(out_dir, "/shared_snp_dt_table_test100.csv")
out_rds <- paste0(out_dir, "/shared_snp_dt_table_test100.rds")
if(!file.exists(out_csv)) file.create(out_csv)
if(!file.exists(out_rds)) file.create(out_rds)


# load melanogaster gds file
mel_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.gds"
sim_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds"
mel_gds <- seqOpen(mel_file)
sim_gds <- seqOpen(sim_file)

# genofile <- sim_file

# genofile <- seqOpen(genofile)
# genofile

# # load simulans gds file
# sim_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds"
# sim_gds <- seqOpen(sim_file)
# sim_gds 

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
get_gds_data <- function(gds, shared, species){
    snp_id <- shared[[paste0("id_", species)]]
    seqSetFilter(gds, variant.id = snp_id)

    message("getting alleles")
    alleles_all <- seqGetData(gds, "allele")
    allele_split <- tstrsplit(alleles_all, ",")

    snp.dt1 <- data.table(
    variant.id = snp_id,
    chr        = shared$chr,
    pos        = shared$pos,
    ref        = allele_split[[1]],
    alt        = allele_split[[2]])

    message("getting annotations")
    ann_all <- seqGetData(gds, "annotation/info/ANN")
    ann_dt <- data.table( variant.id = rep(snp_id, times=ann_all$length), ann = ann_all$data)
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

    # should i be collapsing to one annotation? if so what should i keep?
    # ann_canonical <- ann_dt[, .SD[1], by = variant.id]

    shared_table <- merge(snp.dt1, ann_dt[, .(variant.id, effect, impact, gene, gene_id, feature_type, 
        transcript_id, biotype, in_exon, nt_change, aa_change, aa_pos)], by = "variant.id")
    
    return(shared_table)
}

mel_snp_dt <- build_snp_dt(mel_gds)
sim_snp_dt <- build_snp_dt(sim_gds)

# build shared datatable
shared <- merge(mel_snp_dt, sim_snp_dt, by = c("chr", "pos"), suffixes = c("_mel", "_sim"))
message(nrow(shared), " shared variants")

shared_test <- shared[1:100]
shared_table <- get_gds_data(mel_gds, shared_test, "mel")
# shared_table <- get_gds_data(sim_gds, shared, "sim")
message("saving csv to ", out_csv)
fwrite(shared_table, out_csv)
message("saving rds to ", out_rds)
saveRDS(shared_table, out_rds)

message("complete. ", nrow(shared_table), " variants written.")




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
