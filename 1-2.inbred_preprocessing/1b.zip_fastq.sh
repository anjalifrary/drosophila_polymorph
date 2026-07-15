#!/usr/bin/env bash
#
#SBATCH -J zipFASTQ # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/zip.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/zip.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-184
#ijob -A berglandlab -c10 -p standard --mem=40G

module load gcc htslib


#for pipeline, get sample_dir from array
wd="/scratch/ejy4bu/drosophila/inbred/fastq/PRJNA318623"
# sample_dir="$1"
# sample_dir="${wd}/SRR3585384"
sample_dir=$(find "$wd" -maxdepth 1 -type d -name "SRR*" | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "processing ${sample_dir}"

samp_name=$(basename "$sample_dir")

if [ -f "${sample_dir}/${samp_name}.fastq" ]; then
  gzip ${sample_dir}/${samp_name}.fastq
fi
