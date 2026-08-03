#!/usr/bin/env bash
#
#SBATCH -J be_trim # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/ejy4bu/err_outs/be/trim.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/be/trim.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=0-183%10



### 2. fastp trim single-end mode
module load miniforge && conda activate fastp

if [ ! -f "${SAMPLE_DIR}/fastp/${samp_name}.trimmed.fastq.gz" ]; then
    echo "trimming sample ${samp_name}. "
    fastp \
        --in1 ${fastq} \
        --out1 ${SAMPLE_DIR}/fastp/${samp_name}.trimmed.fastq.gz \
        --json ${SAMPLE_DIR}/fastp/${samp_name}.fastp.json \
    --html ${SAMPLE_DIR}/fastp/${samp_name}.fastp.html \
    --thread 10

    conda deactivate
    conda deactivate
fi