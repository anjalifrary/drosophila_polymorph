#!/usr/bin/env bash
#
#SBATCH -J copy_bam # A single job name for the array
#SBATCH --cpus-per-task=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 
#SBATCH --mem 10G
#SBATCH -o /scratch/ejy4bu/err_outs/be/copy.%A_%a.out  # Std out
#SBATCH -e /scratch/ejy4bu/err_outs/be/copy.%A_%a.out  # Std error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-3666%100

set -euo pipefail

new_dir="/project/berglandlab/alan/be_flies/05.bam"
wd="/scratch/ejy4bu/backyardEvolution/fastq/"
SAMPLE_LIST="/scratch/ejy4bu/backyardEvolution/metadata/allSamples.txt"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
echo "Processing ${sample}"

module load gcc/11.4.0 sratoolkit/3.1.1 
# module load gcc/11.4.0 sratoolkit/3.1.1 aspera-connect/4.2.8

if [ ! -d $new_dir ]; then
  mkdir $new_dir
fi

SAMPLE_DIR="${new_dir}/${sample}"
mkdir -p ${SAMPLE_DIR}

echo "copying $sample..."
for file in \
    "${wd}/${sample}/${sample}.sorted.markdup.clipped.bam" "${SAMPLE_DIR}" \
    "${wd}/${sample}/${sample}.sorted.markdup.clipped.bam.bai" "${SAMPLE_DIR}" \
    "${wd}/${sample}/${sample}.fastp.json" "${SAMPLE_DIR}" \
    "${wd}/${sample}/${sample}.fastp.html" "${SAMPLE_DIR}" \
    "${wd}/${sample}/${sample}.sorted.markdup.metrics.txt" "${SAMPLE_DIR}" \
    "${wd}/${sample}/${sample}.clipOverlapStats.txt" "${SAMPLE_DIR}"
do
    if [ ! -f "$file" ]; then
        echo "ERROR: Missing $file"
        exit 1
    fi
done

echo "copying $sample complete"
