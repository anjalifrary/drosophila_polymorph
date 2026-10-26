#!/usr/bin/env bash
#
#SBATCH -J norm # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/norm_vcf.%A.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/norm_vcf.%A.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

set -euo pipefail
module load gcc/11.4.0
module load bcftools

### dsim Signor
outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/"
ref="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"
in_vcf="${outdir}/dsim3.signor.combined.raw.vcf.gz"
out_vcf="${outdir}/dsim3.signor.combined.norm.vcf.gz"

echo "number records before: "
bcftools index -n $in_vcf


# # ### 6. Normalize vcfs - SIM
# flag -m -both splits all multiallelic records into different rows by type - indel & snp

bcftools norm \
    -f ${ref} \
    -m -both \
    -Oz \
    -o ${out_vcf} \
    ${in_vcf}

echo "normalized"

bcftools index -f ${out_vcf}
echo "indexed. "

echo "number records before: "
bcftools index -n $in_vcf
echo "number records after: "
bcftools index -n $out_vcf

echo "complete."


### 6. Normalize vcfs - MEL

### dmel DGRP2
# outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/"
# ref=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna
# in_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.vcf.gz"
# reheadered_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.vcf.gz"
# prm_chr_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.vcf.gz"
# out_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.vcf.gz"

bcftools reheader \
    -h /scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/DRGP2.dm6.header.txt \
    -o $reheadered_vcf \
    $in_vcf

bcftools index $reheadered_vcf

bcftools view \
    -r 2L,2R,3L,3R,4,X \
    -Oz \
    -o ${prm_chr_vcf} \
    ${reheadered_vcf}

echo "Filtered for primary chr"

bcftools index -f $prm_chr_vcf

bcftools norm \
    -f ${ref} \
    -m -both \
    -Oz \
    -o ${out_vcf} \
    ${prm_chr_vcf}

echo "normalized with -both flag"

bcftools index -f ${out_vcf}
echo "indexed. "

echo "number records before: "
bcftools index -n $in_vcf
echo "number records reheadered: "
bcftools index -n $reheadered_vcf
echo "number records primary chr: "
bcftools index -n $prm_chr_vcf
echo "number records after norm: "
bcftools index -n $out_vcf

echo "complete..."