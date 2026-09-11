library(data.table)
library(stringr)
library(Biostrings)

#######################################
### sim signor

lifted_vcf <- "/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.vcf.gz"

sim_gff <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.gff"

genome_fna <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"
cds_fasta <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim//GCF_016746395.2_Prin_Dsim_3.1_cds_genomic.fna"

sim_dt <- readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.snp_dt_SynMissense.rds")
gds_file <- "/scratch/ejy4bu/drosophila/inbred/gds_files/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.gds"

out_dir <- "/scratch/ejy4bu/drosophila/inbred/snpDT/simMap/"

out_rds <- paste0(out_dir,  "dsim3.signor.all_quality_variants_codonMap.rds")


########

### gff format:
# seqid   source  feature start   end     score   strand  phase   attribute ---   
# sim_2R  Gnomon  CDS     6534278 6534771 .       -       0       ID=cds-XP_016026650.1;Parent=rna-XM_016167191.3;Dbxref=GeneID:6733679,Genbank:XP_016026650.1;Name=XP_016026650.1;gbkey=CDS;gene=LOC6733679;product=dnaJ homolog subfamily C member 13 isoform X1;protein_id=XP_016026650.1

# phase = how many bases to skip at the start of a coding sequence (CDS) feature to reach the first complete codon.

gff <- fread(
     cmd = paste("grep -v '^#'", shQuote(sim_gff)), 
    #  sim_gff, 
     sep="\t", header=F, 
    #  skip = "#",, 
    col.names=c(
        "chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"
    ))

# gff = 391103 rows
# cds rows in gff = 160968

cds_gff <- gff[feature=="CDS", .(
    chr, start, end, strand, phase, attributes
)]

# extract genes and transcripts
cds_gff[, transcript_id := str_extract(attributes,"(?<=Parent=rna-)[^;]+")]
cds_gff[, gene_id := str_extract(attributes, "(?<=Dbxref=GeneID:)[^,;]+")]

### something is weird with transcript ids because sim dt contains transcripts that aren't in cds_gff ?? 
# seems like some transcripts couldn't be mapped to a cds region, even tho they are protein coding... 
transcript_ids <- unique(sim_dt$transcript_id)
cds_keep <- cds_gff[transcript_id%in%transcript_ids]

lift_dt <- fread(
    "/scratch/ejy4bu/drosophila/inbred/snpDT/simMap/liftover_map.tsv",
    col.names = c(
        "chr", "pos",
        "ref", "alt",
        "chr_src", "pos_src",
        "ref_alt_src",
        "flip", "swap"
    )
)

lift_dt[flip == ".", flip := NA_integer_]
lift_dt[swap == ".", swap := NA_integer_]

# fix ref / alt src columns
lift_dt[, c("ref_src", "alt_src") :=
        tstrsplit(ref_alt_src, ",", fixed = TRUE)
]
lift_dt[, ref_alt_src := NULL]

# make integer columns integers
lift_dt[, `:=`(
    pos = as.integer(pos),
    pos_src = as.integer(pos_src),
    flip = as.integer(flip),
    swap = as.integer(swap)
)]
