# library(SeqArray)
# library(data.table)
# library(foreach)
# library(doMC)


# out_dir <- "/scratch/ejy4bu/drosophila/inbred/snpDT/"
# in_rds <- paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt.rds")
# # background_rds <- paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.background_1polymorphicCodon.rds")
# # candidate_rds <- paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.candidate_sharedPolyCodonsOnly.rds")
# out_rds <- paste0(out_dir, "dsim3.signor.DGRP2.source_BCM-HGSC.all_quality_variants_merge_unfilt_codonIDs.rds")


### in bash:
lifted_vcf="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.vcf.gz"
# CHROM = dm6
# POS = dm6
# REF / ALT = dm6
# SRC_CHROM = dsim3
# SRC_POS = dsim3
# SRC_REF_ALT = og dsim3 alleles
# FLIP = if strand was flipped (0 or 1)
# SWAP = if reference / alternate allele ordering changed

# bcftools query \
#   -r 2R:9172292 \
#   -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SRC_CHROM\t%INFO/SRC_POS\t%INFO/SRC_REF_ALT\t%INFO/FLIP\t%INFO/SWAP\t%INFO/EFF\n' \
#   "$liftover_vcf"

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SRC_CHROM\t%INFO/SRC_POS\t%INFO/SRC_REF_ALT\t%INFO/FLIP\t%INFO/SWAP\n' \
    "$lifted_vcf" \
    > /scratch/ejy4bu/anjali/drosophila/inbred/snpDT/simMap/liftover_map.tsv

    
sim_gff="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/ncbi_dataset/data/GCF_016746395.2/GCF_016746395.2_Prin_Dsim_3.1_genomic.gff"

### gff format:
# seqid   source  feature start   end     score   strand  phase   attribute ---   
# sim_2R  Gnomon  CDS     6534278 6534771 .       -       0       ID=cds-XP_016026650.1;Parent=rna-XM_016167191.3;Dbxref=GeneID:6733679,Genbank:XP_016026650.1;Name=XP_016026650.1;gbkey=CDS;gene=LOC6733679;product=dnaJ homolog subfamily C member 13 isoform X1;protein_id=XP_016026650.1


# shared_dt <- readRDS(in_rds)

old_fna="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"

new_fna="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/ncbi_dataset/data/GCF_016746395.2/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna"

cds_fasta="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim//GCF_016746395.2_Prin_Dsim_3.1_cds_genomic.fna"