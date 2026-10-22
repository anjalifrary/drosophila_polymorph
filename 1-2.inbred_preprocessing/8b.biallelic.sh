#!/usr/bin/env bash

#SBATCH -J biallelic               # A single job name for the array
#SBATCH --ntasks-per-node=10    # one core
#SBATCH -N 1                    # on one node
#SBATCH -t 0-10:00              # 10 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/biallelic.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/biallelic.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

set -euo pipefail

module load gcc/11.4.0
module load bcftools

### SIM
# outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/"
# in_vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.vcf.gz"
# out_vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.biallelic.vcf.gz"
# ref="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"


### MEL
ref=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna
indir="/project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/DGRP2/"
outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/"
in_vcf="${indir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.vcf.gz"
out_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.vcf.gz"

echo "Input records:"
bcftools index -n "$in_vcf"
echo "$in_vcf"

bcftools norm \
    --threads 10 \
    -f ${ref} \
    -m +both \
    -Oz \
    -o ${out_vcf} \
    ${in_vcf}

bcftools index -f ${out_vcf}

echo "Merged records: "
bcftools index -n "$out_vcf"

bcftools view \
    --threads 10 \
    -m2 -M2 -v snps \
    -Oz \
    -o ${out_vcf%.vcf.gz}.biallelic.vcf.gz \
    ${out_vcf}

bcftools index -f ${out_vcf%.vcf.gz}.biallelic.vcf.gz


echo "Biallelic records:"
bcftools index -n "${out_vcf%.vcf.gz}.biallelic.vcf.gz"