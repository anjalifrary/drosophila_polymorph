#!/usr/bin/env bash
#
#SBATCH -J norm # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/dest/norm_vcf.%A.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/dest/norm_vcf.%A.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


set -euo pipefail
module load gcc/11.4.0
module load bcftools

# ### DEST mel:
# outdir="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/filtering/"
# in_vcf="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.eff.vcf.gz"
# out_vcf="${outdir}/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.vcf.gz"
# ref=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna

# ### DEST sim:
outdir="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/filtering/"
in_vcf="${outdir}/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.primaryChrs.vcf.gz"
out_vcf="${outdir}/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.primaryChrs.norm.vcf.gz"
ref="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"

echo "number records before: "
bcftools index -n $in_vcf

# # ### 4. Normalize vcfs
# splits all multiallelic sites into biallelic rows 
bcftools norm \
    --threads 10 \
    -f ${ref} \
    -m -both \
    -Oz \
    -o ${out_vcf} \
    ${in_vcf}

echo "normalized with -both flag"

bcftools index -f ${out_vcf}
echo "indexed. "

echo "number of records after:"
bcftools index -n ${out_vcf}

echo "complete..."