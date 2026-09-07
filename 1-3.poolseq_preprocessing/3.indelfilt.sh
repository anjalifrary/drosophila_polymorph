#!/usr/bin/env bash
#
#SBATCH -J indel10bp # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/dest/filt_indel.%A.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/dest/filt_indel.%A.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

set -euo pipefail


module load bcftools
# module load bedtools

# ### DEST mel:
outdir="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/filtering/mel/"
in_vcf="${outdir}/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.vcf.gz"
gap_vcf="${outdir}/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.snpGap10.vcf.gz"
snp_vcf="${outdir}/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.snpGap10.snpsOnly.vcf.gz"

# ### DEST sim:
# outdir="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/filtering/sim/"
# in_vcf="${outdir}/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.norm.vcf.gz"
# gap_vcf="${outdir}/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.norm.snpGap10.vcf.gz"
# snp_vcf="${outdir}/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.norm.snpGap10.snpsOnly.vcf.gz"


echo "filtering via SnpGap"
bcftools filter \
    --SnpGap 10 \
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
    -v snps \
    "$gap_vcf" \
    -Oz \
    -o "$snp_vcf"

bcftools index -t "$snp_vcf"
echo "indexed"

echo "number records before: " 
bcftools index -n $in_vcf

echo "number records after removing sites flanking indels: " 
bcftools index -n $gap_vcf

echo "number records after filtering for snps only: " 
bcftools index -n $snp_vcf

echo "complete"

