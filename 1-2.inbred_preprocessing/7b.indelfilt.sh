#!/usr/bin/env bash
#
#SBATCH -J indel10bp # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/filt_indel.%A.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/filt_indel.%A.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

set -euo pipefail


module load bcftools

# # ### Sim signor: 
# outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/"
# in_vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.vcf.gz"
# gap_vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.snpgap10.vcf.gz"
# snp_vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.vcf.gz"


# ### Mel DGRP2:
# outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/"
# in_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.vcf.gz"
# gap_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.vcf.gz"
# snp_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.vcf.gz"

echo "filtering via SnpGap"
bcftools filter \
    --SnpGap 10 \
    -s SnpGap10 \
    --threads 10 \
    $in_vcf \
    -Oz \
    -o $gap_vcf
echo "completed SnpGap filtering"
bgzip -t "$gap_vcf"
echo "passed test"

bcftools index -t "$gap_vcf"

echo "filtering for snps only"
bcftools view \
    --threads 10 \
    -i 'FILTER!="SnpGap10"' \
    -v snps \
    "$gap_vcf" \
    -Oz \
    -o "$snp_vcf"

bcftools index -t "$snp_vcf"
echo "complete"
