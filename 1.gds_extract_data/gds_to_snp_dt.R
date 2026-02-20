library(SeqArray)
library(data.table)

out_file <- "/scratch/ejy4bu/drosophila/gds_analysis/sim_snp_dt_table.csv"
if(!file.exists(out_file)) file.create(out_file)

# load melanogaster gds file
# mel_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.gds"
sim_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds"
genofile <- sim_file

genofile <- seqOpen(genofile)
genofile

# # load simulans gds file
# sim_file <- "/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds"
# sim_gds <- seqOpen(sim_file)
# sim_gds 


snp.dt <- data.table(
    chr=seqGetData(genofile, "chromosome"),
    pos=seqGetData(genofile, "position"),
    nAlleles=seqGetData(genofile, "$num_allele"),
    id=seqGetData(genofile,"variant.id"))

snp.dt <- snp.dt[nAlleles==2] ##subset to sites with 2 alleles
seqResetFilter(genofile)
seqSetFilter(genofile, variant.id = snp.dt$id)  # biallelic only


message("gettting alleles")
alleles_all <- seqGetData(genofile, "allele")
allele_split <- tstrsplit(alleles_all, ",")

snp.dt1 <- data.table(
    variant.id = snp.dt$id,
    chr        = snp.dt$chr,
    pos        = snp.dt$pos,
    ref        = allele_split[[1]],
    alt        = allele_split[[2]])

message("getting annotations")
ann_all <- seqGetData(genofile, "annotation/info/ANN")

ann_dt <- data.table( variant.id = rep(snp.dt$id, times=ann_all$length), 
    ann = ann_all$data)
ann_split <- tstrsplit(ann_dt$ann, "\\|")

#  $ length: int 15
#  $ data  : chr [1:15] "A|upstream_gene_variant|MODIFIER|CG11023|FBgn0031208|transcript|FBtr0300690|protein_coding||c.-1427T>A|||||1276|" 
# "A|upstream_gene_variant|MODIFIER|CG11023|FBgn0031208|transcript|FBtr0300689|protein_coding||c.-1427T>A|||||1276|" 
# "A|upstream_gene_variant|MODIFIER|CG11023|FBgn0031208|transcript|FBtr0330654|protein_coding||c.-1427T>A|||||1276|" 
# "A|downstream_gene_variant|MODIFIER|l(2)gl|FBgn0002121|transcript|FBtr0078170|protein_coding||c.*4962A>T|||||3586|" ...   
#  - attr(*, "class")= chr "SeqVarDataList"

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

# message("classifying as syn/nonsyn")

# message("intron/exon classification")

message("collapse to canonical transcript...")
ann_canonical <- ann_dt[, .SD[1], by = variant.id]

# join with SNP info
snp_table <- merge(snp.dt1, ann_canonical[, .(variant.id, effect, 
    impact, gene, gene_id, feature_type, transcript_id, biotype, in_exon, 
    nt_change, aa_change, aa_pos)],by = "variant.id")

message("saving to ", out_file)
fwrite(snp_table, out_file)
message("complete. ", nrow(snp_table), " variants written.")
