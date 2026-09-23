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

mel_dt <- readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/DGRP2.source_BCM-HGSC.dm6.snp_dt_SynMissense.rds")

mel_gff <- "/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff"


gff <- fread(
     cmd = paste("grep -v '^#'", shQuote(mel_gff)), 
    #  sim_gff, 
     sep="\t", header=F, 
    #  skip = "#",, 
    col.names=c(
        "chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"
    ))

# map chr names to dm6 normal names 
chr_map <- gff[ feature == "region" & grepl("genome=chromosome", attributes),
    .(
        chr_gff = chr,
        chr = sub(".*Name=([^;]+);.*", "\\1", attributes)
    )
]

gff[
     chr_map,
     chr_dm6 := i.chr,
     on = .(chr = chr_gff)
 ]
    
unique(gff[!is.na(chr_dm6), .(chr, chr_dm6)]) # X, 2L, 2R, 3L, 3R, 4, Y

cds_gff <- gff[
    feature == "CDS" & !is.na(chr_dm6), 
    .(chr_dm6, start, end, strand, phase, attributes)
]

# extract genes and transcripts

mrna_map <- gff[
    feature == "mRNA",
    .(
        transcript_id = str_extract(attributes, "(?<=ID=rna-)[^;]+"),
        transcript_id_fbtr = str_extract(attributes, "FBtr[0-9]+")
    )
]
cds_gff[
    mrna_map,
    on = "transcript_id",
    transcript_id_fbtr := i.transcript_id_fbtr
]

# cds_gff[, transcript_id_fbtr := str_extract(attributes, "FBtr[0-9]+")]
cds_gff[, gene_name := str_extract(attributes, "(?<=;gene=)[^;]+")]
cds_gff[, gene_id_fbgn := str_extract(attributes, "FBgn[0-9]+")]
cds_gff[, .N, by = .(chr_dm6, strand)]


# names(cds_gff)
# [1] "chr_dm6"       "start"         "end"           "strand"        "phase"        
# [6] "attributes"    "transcript_id" "gene_name"     "gene_id_fbgn" 
# > names(shared_dt)
#  [1] "chr"                  "pos"                  "variant.id_mel"      
#  [4] "ref_mel"              "alt_mel"              "af_mel"              
#  [7] "maf_mel"              "n_samps_mel"          "effect_mel"          
# [10] "impact_mel"           "functional_class_mel" "codon_change_mel"    
# [13] "aa_change_mel"        "aa_length_mel"        "gene_mel"            
# [16] "biotype_mel"          "gene_coding_mel"      "transcript_id_mel"   
# [19] "exon_rank_mel"        "genotype_mel"         "variant.id_sim"      
# [22] "ref_sim"              "alt_sim"              "af_sim"              
# [25] "maf_sim"              "n_samps_sim"          "chr_src"             
# [28] "pos_src"              "ref_src"              "alt_src"             
# [31] "flip"                 "swap"                 "effect_sim"          
# [34] "impact_sim"           "functional_class_sim" "codon_change_sim"    
# [37] "aa_change_sim"        "aa_length_sim"        "gene_sim"            
# [40] "biotype_sim"          "gene_coding_sim"      "transcript_id_sim"   
# [43] "exon_rank_sim"        "genotype_sim"         "gene_id_fbgn"        
# [46] "aa_ref_mel"           "aa_alt_mel"           "aa_ref_sim"          
# [49] "aa_alt_sim"           "aa_pos_mel"           "aa_pos_sim"          
# [52] "codon_ref_mel"        "codon_alt_mel"        "codon_ref_sim"       
# [55] "codon_alt_sim"        "dist_prev_mel"        "dist_next_mel"       
# [58] "dist_prev_sim"        "dist_next_sim"        "keep_sim"  

# x = mel_cds
# i = shared_dt
mel_cds <- cds_gff[
    shared_dt[!is.na(variant.id_mel)], 
    on = .(
        chr_dm6 = chr,
        transcript_id_fbtr = transcript_id_mel,
        start <= pos, 
        end >= pos
    ), nomatch = 0,
    .(# GFF information
        variant_id_mel = i.variant.id_mel,
        chr = i.chr,
        pos = i.pos,
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



        ### gff info:
        chr_fromGFF = x.chr_dm6,
        cds_start = x.start,
        cds_end = x.end,
        strand = x.strand,
        phase = x.phase,
        attributes = x.attributes,
        transcript_id_fbtr = x.transcript_id_fbtr,
        transcript_id = x.transcript_id,
        gene_name = x.gene_name,
        gene_id_fbgn_fromGFF = x.gene_id_fbgn,
        gene_id_fbgn = i.gene_id_fbgn
    )
]

mel_cds[, nt_pos_mel := as.integer(sub(".*c\\.([0-9]+).*", "\\1", aa_change_mel))]
mel_cds[, pos_in_codon := ((nt_pos_mel - 1) %% 3) + 1]
# confirm that this agrees with the snpEff uppercase annotation:
mel_cds[, codon_start_nt := ((nt_pos_mel - 1) %/% 3) * 3 + 1]
mel_cds[, codon_end_nt   := codon_start_nt + 2]

mel_cds[, snp_pos_in_codon := regexpr("[A-Z]", sub("/.*", "", codon_change_mel))]
nrow(mel_cds[pos_in_codon != snp_pos_in_codon]) # should be 0 . good!

mel_cds[strand == "+",
    `:=`(
        codon_genomic_start = pos - (pos_in_codon - 1),
        codon_genomic_end   = pos + (3 - pos_in_codon)
    )]

mel_cds[strand == "-",
    `:=`(
        codon_genomic_start = pos - (3 - pos_in_codon),
        codon_genomic_end   = pos + (pos_in_codon - 1)
    )]

### check against reference fna


### check against reference
library(Biostrings)
library(data.table)

ref <- readDNAStringSet("/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna")
mel_cds[, codon_ref_annotation := toupper(sub("/.*", "", codon_change_mel))]
mel_cds[, codon_genomic_raw := NA_character_]

ref_chr <- ref[trimws(names(ref)) %in% chr_map$chr]

names(ref_chr) <- trimws(names(ref_chr))
for (chr in unique(mel_cds$chr)) {
    idx <- which(mel_cds$chr == chr)
    seq <- as.character(ref_chr[[chr]])

    mel_cds$codon_genomic_raw[idx] <- substring(
        seq,
        mel_cds$codon_genomic_start[idx],
        mel_cds$codon_genomic_end[idx]
    )
}


mel_cds[strand == "+", codon_from_genome := codon_genomic_raw]
mel_cds[strand == "-", codon_from_genome := as.character( reverseComplement(DNAStringSet(as.character(codon_genomic_raw))))]

nrow(mel_cds[codon_from_genome == codon_ref_annotation])
nrow(mel_cds[codon_from_genome != codon_ref_annotation])

table(mel_cds$codon_from_genome == mel_cds$codon_ref_annotation)

# codon_genomic_start = genomic coord of 1st base on reference strand
#              _end = 3rd base 
# codon_genomic_raw = the actual codon (eg: CCT) read left to right on reference strand (+)
# codon_from_genome = genomic codon converted into 5' to 3' (reverse complmented for - strand)
# codon_ref_annotation = reference codon annotated by snpeff (from codon_change_sim)


### verified that codon calculated from reference matches that calculated from snp table



### next: create a codon id of sim-based coordinates: chr:[start]:[end]
mel_cds[, codon_id := paste(
    chr,
    codon_genomic_start,
    codon_genomic_end,
    sep = ":"
)]

mel_cds[, n_snps_in_codon := .N, by = codon_id]

mel_1snp_perCodon <- mel_cds[n_snps_in_codon ==1]

saveRDS(mel_cds, "/scratch/ejy4bu/drosophila/inbred/snpDT/DGRP2.source_BCM-HGSC.dm6.snp_dt_SynMissense_update09_17_26.rds")
saveRDS(mel_1snp_perCodon, "/scratch/ejy4bu/drosophila/inbred/snpDT/DGRP2.source_BCM-HGSC.dm6.snp_dt_SynMissense_update09_17_26_1snp_perCodon.rds")

shared_dt[, keep_mel := as.character(variant.id_mel %in% mel_1snp_perCodon$variant_id_mel)]

shared_dt[, keep_mel := fifelse(
    is.na(variant.id_mel),
    NA_character_,
    as.character(variant.id_mel %in% mel_1snp_perCodon$variant_id_mel)
)]

# this annotates all rows TRUE or FALSE depending on if they are in the sim dt 
    # so FALSE includes rows that didn't pass multiple snp per codon filter AND 
    #   those that have no sim annotation at that row
# the number of TRUE rows == number of rows in sim_1snp_perCodon



# there are 6197 rows that are missing from mel cds (annotated as FALSE but have variant id mel)
missing_from_mel_cds <- shared_dt[
    keep_mel == FALSE & !is.na(variant.id_mel) &
    !(variant.id_mel %in% mel_cds$variant_id_mel)
]

nrow(missing_from_mel_cds)
# they are all unassigned transcripts.. 

shared_dt[variant.id_mel %in% missing_from_mel_cds$variant.id_mel, keep_mel := "UNASSIGNED"]
table(shared_dt$keep_mel, useNA="ifany")

saveRDS(shared_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.rds")



############################################################################################
### check if any of the sim liftover results show up as in the wrong spot for mel 
    # in order to determine same / different codon snps

shared_dt <- readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.rds")
mel_1snp_perCodon <- readRDS("/scratch/ejy4bu/drosophila/inbred/snpDT/DGRP2.source_BCM-HGSC.dm6.snp_dt_SynMissense_update09_17_26_1snp_perCodon.rds")
sim_1snp_perCodon <- readRDS(paste0(out_dir, "sim_1snp_perCodon.rds"))

### add dm6 codon id for sim variants:
sim_1snp_perCodon[
    strand == "+",
    `:=`(
        dm6_codon_start = pos - (snp_pos_in_codon - 1),
        dm6_codon_end   = pos + (3 - snp_pos_in_codon)
    )
]

sim_1snp_perCodon[
    strand == "-",
    `:=`(
        dm6_codon_start = pos - (3 - snp_pos_in_codon),
        dm6_codon_end   = pos + (snp_pos_in_codon - 1)
    )
]

sim_1snp_perCodon[
    ,
    dm6_codon_id := paste(
        chr,
        dm6_codon_start,
        dm6_codon_end,
        sep = ":"
    )
]

### rename mel codon id for dm6 coordinates
mel_1snp_perCodon[
    ,
    dm6_codon_id := codon_id
]

### ID codons where 2 sim dsim3 codons map to 1 dm6 codon
dup_dm6_codons <- sim_1snp_perCodon[
    ,
    .N,
    by = dm6_codon_id
][N > 1]


shared_dt[
    variant.id_sim %in% sim_1snp_perCodon[
        dm6_codon_id %in% dup_dm6_codons$dm6_codon_id,
        variant_id_sim
    ],
    keep_sim := "MAP2TO1"
]

table(shared_dt$keep_sim, useNA = "ifany")
shared_dt_all <- shared_dt[
    sim_1snp_perCodon,
    on = .(variant.id_sim = variant_id_sim),
    `:=`(
        strand_sim = i.strand,
        snp_pos_in_codon_sim = i.snp_pos_in_codon,
        codon_id_dsim3 = i.codon_id,
        codon_id_dm6_sim = i.dm6_codon_id
    )
]

shared_dt_all <- shared_dt_all[
    mel_1snp_perCodon,
    on = .(variant.id_mel = variant_id_mel),
    `:=`(
        strand_mel = i.strand,
        snp_pos_in_codon_mel = i.snp_pos_in_codon,
        codon_id_dm6_mel = i.dm6_codon_id,
        gene_id_fbgn_fromGFF=i.gene_id_fbgn_fromGFF, 
        gene_id_fbgn=i.gene_id_fbgn
    )
]
    
# (, .(strand, snp_pos_in_codon, codon_id, dm6_codon_id),)]
shared_dt_all[, classification := NA_character_] # add empty classification colun

shared_dt_cleaner <- shared_dt_all[, .(
    ### meta stuff 
    classification,  
    variant.id_mel, variant.id_sim,
    chr_dm6=chr, pos_dm6=pos,
    chr_dsim3=chr_src, pos_dsim3=pos_src,
    
    ### nt stuff:
    ref_mel, alt_mel, 
    ref_sim_dm6=ref_sim, alt_sim_dm6=alt_sim, 
    ref_sim_dsim3=ref_src, alt_sim_dsim3=alt_src, 

    ### codon stuff:
    codon_id_dm6_mel, codon_id_dm6_sim, codon_id_dsim3,
    codon_ref_mel, codon_alt_mel, codon_ref_sim, codon_alt_sim,
    snp_pos_in_codon_mel, snp_pos_in_codon_sim,

    ### amino acid stuff:
    aa_pos_mel, aa_pos_sim,
    aa_ref_mel, aa_alt_mel, aa_ref_sim, aa_alt_sim,

    ### gene / transcript stuff:
    gene_mel, gene_sim, 
    gene_id_fbgn,
    gene_id_fbgn_fromGFF,

    transcript_id_mel, transcript_id_sim,

    ### frequency stuff:
    n_samps_mel, n_samps_sim, 
    af_mel, maf_mel, af_sim, maf_sim, 
    effect_mel, effect_sim, 

    ### extras:
    flip, swap, 
    strand_mel, strand_sim,
    keep_mel, keep_sim
)]

saveRDS(shared_dt_cleaner, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.cleaner.rds")

### clean up rows where keep_mel / sim != TRUE

filtered_dt <- shared_dt_cleaner
filtered_dt[keep_mel != "TRUE", variant.id_mel := NA_integer_]
filtered_dt[keep_sim != "TRUE", variant.id_sim := NA_integer_]

filtered_dt <- filtered_dt[!is.na(variant.id_mel) | !is.na(variant.id_sim)]
### mel cols:
mel_cols <- c(
    "ref_mel", "alt_mel", 
    "codon_id_dm6_mel", 
    "codon_ref_mel", "codon_alt_mel",
    "snp_pos_in_codon_mel", 
    "aa_pos_mel", 
    "aa_ref_mel", "aa_alt_mel", 
    "gene_mel", "gene_id_fbgn_fromGFF", 
    "transcript_id_mel", 
    "n_samps_mel", 
    "af_mel", "maf_mel", 
    "effect_mel", 
    "strand_mel"
)
### sim cols:
sim_cols <- c(
    "chr_dsim3", "pos_dsim3", 
    "ref_sim_dm6", "alt_sim_dm6", "ref_sim_dsim3", "alt_sim_dsim3",
    "codon_id_dm6_sim", "codon_id_dsim3",
    "codon_ref_sim", "codon_alt_sim",
    "snp_pos_in_codon_sim", 
    "aa_pos_sim",
    "aa_ref_sim", "aa_alt_sim",
    "gene_sim",
    "transcript_id_sim",
    "n_samps_sim", 
    "af_sim", "maf_sim",
    "effect_sim", 
    "flip", "swap", 
    "strand_sim"
)

filtered_dt[
    is.na(variant.id_mel),
    (mel_cols) := lapply(.SD, function(x) {
        if (is.character(x)) NA_character_ else NA
    }),
    .SDcols = mel_cols
]

filtered_dt[
    is.na(variant.id_sim),
    (sim_cols) := lapply(.SD, function(x) {
        if (is.character(x)) NA_character_ else NA
    }),
    .SDcols = sim_cols
]

saveRDS(filtered_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.cleaner.1snpCodon.rds")


### working off filtered_dt


### checking if there are codons that overlap (different reading frame)
mel_codons <- unique(
    filtered_dt[
        !is.na(variant.id_mel),
        .(
            codon_id = codon_id_dm6_mel,
            chr = chr_dm6,
            start = as.integer(tstrsplit(codon_id_dm6_mel, ":", fixed = TRUE)[[2]]),
            end = as.integer(tstrsplit(codon_id_dm6_mel, ":", fixed = TRUE)[[3]])
        )
    ]
)

sim_codons <- unique(
    filtered_dt[
        !is.na(variant.id_sim),
        .(
            codon_id = codon_id_dm6_sim,
            chr = chr_dm6,
            start = as.integer(tstrsplit(codon_id_dm6_sim, ":", fixed = TRUE)[[2]]),
            end = as.integer(tstrsplit(codon_id_dm6_sim, ":", fixed = TRUE)[[3]])
        )
    ]
)

overlap <- mel_codons[
    sim_codons,
    on = .(
        chr,
        start <= end,
        end >= start
    ),
    nomatch = 0,
    .(
        codon_id_mel=x.codon_id,
        codon_id_sim=i.codon_id
    ),
    allow.cartesian = TRUE
]

shifted_overlap <- overlap[
    codon_id_mel != codon_id_sim
]


filtered_dt[
    codon_id_dm6_mel %in% shifted_overlap$codon_id_mel,
    keep_mel := "OVERLAP"
]

filtered_dt[
    codon_id_dm6_sim %in% shifted_overlap$codon_id_sim,
    keep_sim := "OVERLAP"
]



### inspecting # snps per codons

mel_pos <- filtered_dt[
    !is.na(variant.id_mel),
    .(codon_id = codon_id_dm6_mel,
      mel_pos = snp_pos_in_codon_mel)
]

sim_pos <- filtered_dt[
    !is.na(variant.id_sim),
    .(codon_id = codon_id_dm6_sim,
      sim_pos = snp_pos_in_codon_sim)
]

shared <- merge(mel_pos, sim_pos, by = "codon_id")

shared[, same_position := mel_pos == sim_pos]

table(shared$same_position,useNA="ifany")

shared[, tsp_samePos := same_position]

tsp_map <- unique(shared[, .(
    codon_id = codon_id,
    tsp_samePos
)])

filtered_dt[
    tsp_map,
    on = .(
        codon_id_dm6_mel = codon_id
    ),
    tsp_samePos := i.tsp_samePos
]

filtered_dt[
    tsp_map,
    on = .(
        codon_id_dm6_sim = codon_id
    ),
    tsp_samePos := i.tsp_samePos
]

table(filtered_dt$tsp_samePos, useNA="ifany") 
# TRUE = number of variants with same position variant
# FALSE = 2x number of codons with diff position snps
# NA = species specific variants OR codon id did not match (reading frame was off)


saveRDS(filtered_dt, "/scratch/ejy4bu/drosophila/inbred/snpDT/dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.annotatedSim.annotatedMel.cleaner.1snpCodon.rds")
