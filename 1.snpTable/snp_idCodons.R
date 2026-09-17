library(data.table)
library(stringr)
library(Biostrings)

#######################################
### sim signor

lifted_vcf <- "/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.vcf.gz"

sim_gff <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.gff"

genome_fna <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"
cds_fasta <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_cds_genomic.fna"


sim_dt <- readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.snp_dt_SynMissense.rds")
gds_file <- "/scratch/ejy4bu/drosophila/inbred/gds_files/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.gds"
shared_dt <-readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.rds")

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

cds_gff[, sim_chr := chr]
cds_gff[, chr := sub("^sim_", "", sim_chr)]

sim_cds <- cds_gff[
    shared_dt[!is.na(variant.id_sim)], 
    on = .(
        sim_chr = chr_src,
        transcript_id = transcript_id_sim,
        start <= pos_src, 
        end >= pos_src
    ), nomatch = 0,
    .(# GFF information
        sim_chr = x.sim_chr,
        cds_start = x.start,
        cds_end = x.end,
        strand = x.strand,
        phase = x.phase,
        attributes = x.attributes,
        transcript_id = x.transcript_id,
        gene_id = x.gene_id,

        # shared_dt information
        chr = i.chr,
        pos = i.pos,

        variant_id_mel = i.variant.id_mel,
        ref_mel = i.ref_mel,
        alt_mel = i.alt_mel,
        af_mel = i.af_mel,
        maf_mel = i.maf_mel,
        effect_mel = i.effect_mel,
        codon_change_mel = i.codon_change_mel,
        aa_change_mel = i.aa_change_mel,
        aa_pos_mel = i.aa_pos_mel,

        variant_id_sim = i.variant.id_sim,
        ref_sim = i.ref_sim,
        alt_sim = i.alt_sim,
        af_sim = i.af_sim,
        maf_sim = i.maf_sim,
        effect_sim = i.effect_sim,
        codon_change_sim = i.codon_change_sim,
        aa_change_sim = i.aa_change_sim,
        aa_pos_sim = i.aa_pos_sim,

        # simulans coordinate/liftover information
        transcript_id_sim = i.transcript_id_sim,
        chr_src = i.chr_src,
        pos_src = i.pos_src,
        ref_src = i.ref_src,
        alt_src = i.alt_src,
        flip = i.flip,
        swap = i.swap
    )
]

sim_cds[, nt_pos_sim := as.integer(sub(".*c\\.([0-9]+).*", "\\1", aa_change_sim))]
sim_cds[, pos_in_codon := ((nt_pos_sim - 1) %% 3) + 1]
# confirm that this agrees with the snpEff uppercase annotation:
sim_cds[, codon_start_nt := ((nt_pos_sim - 1) %/% 3) * 3 + 1]
sim_cds[, codon_end_nt   := codon_start_nt + 2]

sim_cds[, snp_pos_in_codon := regexpr("[A-Z]", sub("/.*", "", codon_change_sim))]
nrow(sim_cds[pos_in_codon != snp_pos_in_codon]) # should be 0 . good!

sim_cds[strand == "+",
    `:=`(
        codon_genomic_start = pos_src - (pos_in_codon - 1),
        codon_genomic_end   = pos_src + (3 - pos_in_codon)
    )]

sim_cds[strand == "-",
    `:=`(
        codon_genomic_start = pos_src - (3 - pos_in_codon),
        codon_genomic_end   = pos_src + (pos_in_codon - 1)
    )]

### check against reference
library(Biostrings)
library(data.table)

ref <- readDNAStringSet("/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna")
sim_cds[, codon_ref_annotation := toupper(sub("/.*", "", codon_change_sim))]

for (chr in unique(sim_cds$sim_chr)) {
    
    idx <- which(sim_cds$sim_chr == chr)
    
    seq <- as.character(ref[[chr]])
    
    sim_cds$codon_genomic_raw[idx] <- substring(
        seq,
        sim_cds$codon_genomic_start[idx],
        sim_cds$codon_genomic_end[idx]
    )
}


sim_cds[strand == "+", codon_from_genome := codon_genomic_raw]
sim_cds[strand == "-", codon_from_genome := as.character( reverseComplement(DNAStringSet(as.character(codon_genomic_raw))))]

nrow(sim_cds[codon_from_genome == codon_ref_annotation])
nrow(sim_cds[codon_from_genome != codon_ref_annotation])

table(sim_cds$codon_from_genome == sim_cds$codon_ref_annotation)

# codon_genomic_start = genomic coord of 1st base on reference strand
#              _end = 3rd base 
# codon_genomic_raw = the actual codon (eg: CCT) read left to right on reference strand (+)
# codon_from_genome = genomic codon converted into 5' to 3' (reverse complmented for - strand)
# codon_ref_annotation = reference codon annotated by snpeff (from codon_change_sim)


### verified that codon calculated from reference matches that calculated from snp table

### next: create a codon id of sim-based coordinates: chr:[start]:[end]
sim_cds[, codon_id := paste(
    sim_chr,
    codon_genomic_start,
    codon_genomic_end,
    sep = ":"
)]

sim_cds[, n_snps_in_codon := .N, by = codon_id]

sim_1snp_perCodon <- sim_cds[n_snps_in_codon ==1]

saveRDS(sim_cds, out_rds)
saveRDS(sim_1snp_perCodon, paste0(out_dir, "sim_1snp_perCodon.rds"))

shared_dt[, keep_sim := as.character(variant.id_sim %in% sim_1snp_perCodon$variant_id_sim)]

shared_dt[, keep_sim := fifelse(
    is.na(variant.id_sim),
    NA_character_,
    as.character(variant.id_sim %in% sim_1snp_perCodon$variant_id_sim)
)]

# this annotates all rows TRUE or FALSE depending on if they are in the sim dt 
    # so FALSE includes rows that didn't pass multiple snp per codon filter AND 
    #   those that have no sim annotation at that row
# the number of TRUE rows == number of rows in sim_1snp_perCodon



# there are 88 rows that are missing from sim cds (annotated as FALSE but have variant id sim)
missing_from_sim_cds <- shared_dt[
    keep_sim == FALSE & !is.na(variant.id_sim) &
    !(variant.id_sim %in% sim_cds$variant_id_sim)
]

nrow(missing_from_sim_cds)
# they are all unassigned transcripts.. 

shared_dt[variant.id_sim %in% missing_from_sim_cds$variant.id_sim, keep_sim := "UNASSIGNED"]
table(shared_dt$keep_sim, useNA="ifany")

saveRDS(shared_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.rds")


############################################################################################
### repeat for mel: 1 snp per codon




############################################################################################
### check if any of the sim liftover results show up as in the wrong spot for mel 
    # in order to determine same / different codon snps





### TRASH:

sim_cds <- cds_gff[sim_dt, on=.(
    sim_chr = chr_src,
    transcript_id = transcript_id,
    start <= pos,
    end >= pos
), nomatch = 0, 
.(
        sim_chr = x.sim_chr,
        cds_start = x.start,
        cds_end = x.end,
        strand = x.strand,
        phase = x.phase,
        attributes = x.attributes,
        transcript_id = x.transcript_id,
        gene_id = x.gene_id,

        # SNP information
        snp_pos = i.pos,
        ref = i.ref,
        alt = i.alt,
        af = i.af,
        maf = i.maf,
        n_samps = i.n_samps,
        pos_src = i.pos_src,
        ref_src = i.ref_src,
        alt_src = i.alt_src,
        flip = i.flip,
        swap = i.swap,
        effect = i.effect,
        impact = i.impact,
        functional_class = i.functional_class,
        codon_change = i.codon_change,
        aa_change = i.aa_change,
        aa_length = i.aa_length,
        gene = i.gene,
        biotype = i.biotype,
        gene_coding = i.gene_coding,
        exon_rank = i.exon_rank,
        genotype = i.genotype
    )
]

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
