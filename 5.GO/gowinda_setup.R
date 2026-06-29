library(data.table)

# requires four files:
    # background snps in .txt file of chr | pos
    # candidate snps - subset
    # GTF annotation file - make sure I'm using the right one
    # GO gene sets

dir <- "/scratch/ejy4bu/drosophila/gowinda/MAF5/"

### background snps
rds <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_dt_analysis/merged_tables/quality/all_quality_variants_merge_unfilt.rds")
total_snp <- unique(rds[, .(chr, pos)])
fwrite(total_snp, paste0(dir, "background_all_snps.txt"), sep="\t", col.names=FALSE)

### filtered background snps
voi <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/subset_qualVar_ofInterest_MAF5_06-18-2026.rds")
candidate_snp <- unique(voi[, .(chr, pos)])
fwrite(candidate_snp, paste0(dir, "background_classed_snps.txt"), sep="\t", col.names=FALSE)

### candidate snps by classes
# fwrite(unique(voi[classification=="A", .(chr, pos)]), paste0(dir, "candidate_snp_A.txt"), sep="\t", col.names=FALSE)
# fwrite(unique(voi[classification=="B", .(chr, pos)]), paste0(dir, "candidate_snp_B.txt"), sep="\t", col.names=FALSE)
fwrite(unique(voi[classification%in%c("A", "B"), .(chr, pos)]), paste0(dir, "candidate_snp_AB.txt"), sep="\t", col.names=FALSE)
# fwrite(unique(voi[classification%in%c("F", "G", "O", "P"), .(chr, pos)]), paste0(dir, "candidate_snp_FGOP.txt"), sep="\t", col.names=FALSE)
# fwrite(unique(voi[classification%in%c("X", "Y"), .(chr, pos)]), paste0(dir, "candidate_snp_XY.txt"), sep="\t", col.names=FALSE)
fwrite(unique(voi[classification%in%c("F", "G", "O", "P", "X", "Y"), .(chr, pos)]), paste0(dir, "candidate_snp_FGOPXY.txt"), sep="\t", col.names=FALSE)
fwrite(unique(voi[classification%in%c("A", "B", "F", "G", "O", "P", "X", "Y"), .(chr, pos)]), paste0(dir, "candidate_snp_ABFGOPXY.txt"), sep="\t", col.names=FALSE)

### gtf annotation file
gtf_file <- fread("/scratch/ejy4bu/drosophila/gowinda/dmel-all-r6.67.gtf")

### GO file
gaf <- fread(
    "/scratch/ejy4bu/drosophila/gowinda/gene_association.fb",
    sep="\t", header=FALSE, fill=TRUE, skip=5
)

# Remove comment lines
# gaf <- gaf[!grepl("^!", V1)]
gowinda_go <- gaf[grepl("^FBgn", V2) & grepl("^GO:", V5),
    .(
        TERM  = unique(V5)[1],
        genes = paste(unique(V2), collapse=" ")
    ), by=V5]
    
fwrite(gowinda_go, "/scratch/ejy4bu/drosophila/gowinda/flybase_gaf_go.txt",
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

fwrite(gowinda_go, "/scratch/ejy4bu/drosophila/gowinda/flybase_go.txt",
    sep = "\t", col.names = FALSE, quote = FALSE)

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