#!/usr/bin/env bash
#
#SBATCH -J norm # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/norm_vcf.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/norm_vcf.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

set -euo pipefail

module load bcftools
module load gcc/11.4.0

outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/"
# ref="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"
ref=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna

vcf="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2.source_BCM-HGSC.dm6.final.vcf.gz"

echo "number rows before: "
bcftools index -n $vcf

bcftools index $vcf
bcftools view \
    -r 2L,2R,3L,3R,4,X \
    -Oz \
    -o "/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2.source_BCM-HGSC.dm6.final.primaryChr.vcf.gz" \
    ${vcf}

echo "Filtered for primary chr"

# ### 4. Normalize vcfs
# # what should -m flag be?? +both merges separate rows into 1 multiallelic row by type (snp vs indel)
# bcftools norm \
#     -f ${ref} \
#     -m +both \
#     -Oz \
#     -o ${outdir}/dsim3.signor.combined.norm.vcf.gz \
#     ${outdir}/dsim3.signor.combined.raw.vcf.gz

# # filters for biallelic snps (remove indels, multiallelic snps)
# bcftools view \
#     -v snps \
#     -m2 -M2 \
#     -Oz \
#     -o ${outdir}/dsim3.signor.combined.norm.biallelic.snpsOnly.vcf.gz \
#     ${outdir}/dsim3.signor.combined.norm.vcf.gz

# bcftools index ${outdir}/dsim3.signor.combined.norm.biallelic.snpsOnly.vcf.gz

# mel DGRP vcf:

### 4. Normalize vcfs
# what should -m flag be?? +both merges separate rows into 1 multiallelic row by type (snp vs indel)
bcftools norm \
    -f ${ref} \
    -m +both \
    -Oz \
    -o ${outdir}/DGRP2.source_BCM-HGSC.dm6.final.norm.vcf.gz \
    ${outdir}/DGRP2.source_BCM-HGSC.dm6.final.primaryChr.vcf.gz

echo "normalized"

# filters for biallelic snps (remove indels, multiallelic snps)
bcftools view \
    -v snps \
    -m2 -M2 \
    -Oz \
    -o ${outdir}/DGRP2.source_BCM-HGSC.dm6.final.norm.biallelic.snpsOnly.vcf.gz \
    ${outdir}/DGRP2.source_BCM-HGSC.dm6.final.norm.vcf.gz

echo "filtered for biallelic snps"

echo "number rows after: "
bcftools index -n ${outdir}/DGRP2.source_BCM-HGSC.dm6.final.norm.biallelic.snpsOnly.vcf.gz

bcftools index ${outdir}/DGRP2.source_BCM-HGSC.dm6.final.norm.biallelic.snpsOnly.vcf.gz

echo "indexed. complete..."