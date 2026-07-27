#!/usr/bin/env bash
#
#SBATCH -J genotype # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/genotype.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/genotype.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=0-5


set -euo pipefail

module load gcc/11.4.0
module load gatk
module load bcftools
module load htslib

ref="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"
outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/"
mkdir -p ${outdir}
gvcf_dir="/scratch/ejy4bu/drosophila/inbred/fastq/PRJNA318623"

#  /project/berglandlab/alan/privatePolymorphisms/simulans/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.vcf.gz | \
chroms=("sim_2L" "sim_2R" "sim_3L" "sim_3R" "sim_4" "sim_X")
# chroms=("2L" "2R" "3L" "3R" "4" "X")
chr=${chroms[$SLURM_ARRAY_TASK_ID]}

JAVAMEM=80G
CPU=10

echo "Processing Chromosome: ${chr}"

tmp=/scratch/ejy4bu/tmp/gatk
mkdir -p "$tmp"

# fix for separate chr outputs from genomicsdbimport:
### Joint Genotyping
gatk --java-options "-Xmx${JAVAMEM} -Djava.io.tmpdir=$tmp" GenotypeGVCFs \
    -R ${ref} \
    -V gendb://${outdir}/dsim_genomicsdb_${chr} \
    -O ${outdir}/dsim3.signor.${chr}.raw.vcf.gz


# # # merge all vcfs: (do outside of array job)
# bcftools concat \
#     -Oz \
#     -o ${outdir}/dsim3.signor.combined.raw.vcf.gz \
#     ${outdir}/dsim3.signor.sim_2L.raw.vcf.gz \
#     ${outdir}/dsim3.signor.sim_2R.raw.vcf.gz \
#     ${outdir}/dsim3.signor.sim_3L.raw.vcf.gz \
#     ${outdir}/dsim3.signor.sim_3R.raw.vcf.gz \
#     ${outdir}/dsim3.signor.sim_4.raw.vcf.gz \
#     ${outdir}/dsim3.signor.sim_X.raw.vcf.gz

# bcftools index ${outdir}/dsim3.signor.combined.raw.vcf.gz