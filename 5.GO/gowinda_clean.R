
### Gowinda wiki: https://sourceforge.net/p/gowinda/wiki/Manual/
### FOR LAUREN

library(SeqArray)
library(data.table)
library(doMC)
registerDoMC(16)

# requires four files:
    # background snps in .txt file format: chr pos (tab-delim)
    # candidate snps - subset
    # GTF annotation file       "/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/dmel-all-r6.67.gtf"
    # GO gene sets              "/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/flybase_gaf_go.txt"

### GTF annotation file
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


### making background and candidate snp table
fwrite(candidate_snp_dt[, .(chr, pos)], 
    paste0(paste0(candidate_file, ".txt")), 
    sep="\t", col.names=FALSE)

fwrite(background_snp_dt[, .(chr, pos)], 
    paste0(paste0(background_file, ".txt")), 
    sep="\t", col.names=FALSE)


###############################################################
### BASH SCRIPT TO RUN


# gtf_file=/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/dmel-all-r6.67.gtf
# go_file=/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/flybase_gaf_go.txt
# jar_file=/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/Gowinda-1.12.jar


# java -Xmx8g -jar $jar_file \
#   --snp-file $background \
#   --candidate-snp-file $candidate \
#   --gene-set-file $go_file \
#   --annotation-file $gtf_file \
#   --simulations 1000000 \
#   --gene-definition gene \
#   --threads 10 \
#   --mode gene \
#   --output-file $output

# gene mode is more conservative 
# snp mode is less conservative, more enrichment