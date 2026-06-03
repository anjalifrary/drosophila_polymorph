

library(data.table)

# requires four files:
    # background snps in .txt file of chr | pos
    # candidate snps - subset
    # GTF annotation file - make sure I'm using the right one
    # GO gene sets

mel_rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/species_rdsFiles/mel_filtered_eff_snp_dt.rds")
total_snp <- unique(mel_rds[, .(chr, pos)])
fwrite(total_snp, "/scratch/ejy4bu/drosophila/gowinda/total_snp.txt", sep="\t", col.names=FALSE)

voi <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_classed_geva.rds")
candidate_snp <- unique(voi[, .(chr, pos)])
fwrite(candidate_snp, "/scratch/ejy4bu/drosophila/gowinda/candidate_snp.txt", sep="\t", col.names=FALSE)

gtf_file <- fread("/scratch/ejy4bu/drosophila/gowinda/dmel-all-r6.67.gtf")

gaf <- fread("/scratch/ejy4bu/drosophila/gowinda/gene_association.fb",
    sep = "\t", header = FALSE)

# gowinda_go <- unique(gaf[, .(GO = V5, Gene = V2)])
gowinda_go <- gaf[, .(genes = paste(unique(V2), collapse=";")), by=V5]


fwrite(gowinda_go,"/scratch/ejy4bu/drosophila/gowinda/flybase_go.txt",
    sep = "\t", col.names = FALSE)

fwrite(unique(voi[classification=="B", .(chr, pos)]), "/scratch/ejy4bu/drosophila/gowinda/candidate_snp_B.txt", sep="\t", col.names=FALSE)

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
    genes = paste(unique(FLYBASE), collapse = ";")
), by = GO]

fwrite(gowinda_go, "/scratch/ejy4bu/drosophila/gowinda/flybase_go.txt",
    sep = "\t", col.names = FALSE)

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