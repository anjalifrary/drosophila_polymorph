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


module load gcc/11.4.0
module load bcftools
module load gatk
module load htslib

# ### mel:
# ref=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna
# outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/"
# in_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.norm.vcf.gz"
# out_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.norm.filtered.vcf.gz"

### sim
ref="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"
outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/"
in_vcf="${outdir}/dsim3.signor.combined.norm.vcf.gz"
out_vcf="${outdir}/dsim3.signor.combined.norm.filtered.vcf.gz"

### 5. Filter variants (used GATK best practices hard filters.. ?)
echo "Filtering..."
gatk VariantFiltration \
    -R ${ref} \
    -V ${in_vcf} \
    -O ${out_vcf} \
    --filter-name "QD2" \
    --filter-expression "QD < 2.0" \
    --filter-name "FS60" \
    --filter-expression "FS > 60.0" \
    --filter-name "MQ40" \
    --filter-expression "MQ < 40.0" \
    --filter-name "MQRankSum125" \
    --filter-expression "MQRankSum < -12.5" \
    --filter-name "ReadPosRankSum8" \
    --filter-expression "ReadPosRankSum < -8.0"
