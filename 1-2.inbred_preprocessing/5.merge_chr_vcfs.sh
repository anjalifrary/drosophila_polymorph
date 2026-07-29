#!/usr/bin/env bash
#
#SBATCH -J merge # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/merge.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/merge.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

set -euo pipefail

module load bcftools

outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/"


# # merge all vcfs: (do outside of array job)
bcftools concat \
    -Oz \
    -o ${outdir}/dsim3.signor.combined.raw.vcf.gz \
    ${outdir}/dsim3.signor.sim_2L.raw.vcf.gz \
    ${outdir}/dsim3.signor.sim_2R.raw.vcf.gz \
    ${outdir}/dsim3.signor.sim_3L.raw.vcf.gz \
    ${outdir}/dsim3.signor.sim_3R.raw.vcf.gz \
    ${outdir}/dsim3.signor.sim_4.raw.vcf.gz \
    ${outdir}/dsim3.signor.sim_X.raw.vcf.gz

bcftools index ${outdir}/dsim3.signor.combined.raw.vcf.gz