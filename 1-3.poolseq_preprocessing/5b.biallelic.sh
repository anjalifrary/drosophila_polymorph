#!/usr/bin/env bash

#SBATCH -J norm               # A single job name for the array
#SBATCH --ntasks-per-node=10    # one core
#SBATCH -N 1                    # on one node
#SBATCH -t 0-10:00              # 10 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/norm.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/norm.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

set -euo pipefail

module load bcftools
module load htslib
module load gcc/11.4.0


## mel:
# outdir="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/filtering/mel/"
# in_vcf="${outdir}/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.snpGap10.snpsOnly.repeatmasked.wmdust.ann.eff.vcf.gz"
# out_vcf="${outdir}/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.snpGap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.vcf.gz"
# ref=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna


### sim:
outdir="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/filtering/sim/"
out_vcf="${outdir}/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.norm.snpGap10.snpsOnly.repeatmasked.wmdust.vcf.gz"
out_vcf="${outdir}/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.norm.snpGap10.snpsOnly.repeatmasked.wmdust.merged.vcf.gz"
ref="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"


echo "Input records:"
bcftools index -n "$in_vcf"

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
    -m2 -M2 -v snps \
    -Oz \
    -o ${out_vcf%.vcf.gz}.biallelic.vcf.gz \
    ${out_vcf}

bcftools index -f ${out_vcf%.vcf.gz}.biallelic.vcf.gz


echo "Biallelic records:"
bcftools index -n "${out_vcf%.vcf.gz}.biallelic.vcf.gz"