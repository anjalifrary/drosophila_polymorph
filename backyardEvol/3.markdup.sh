#!/usr/bin/env bash
#
#SBATCH -J be_dup # A single job name for the array
#SBATCH --cpus-per-task=4
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 30G
#SBATCH -o /scratch/ejy4bu/err_outs/be/markdup.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/be/markdup.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


set -euo pipefail

### to Run
# sbatch --array=1-$(wc -l < /scratch/ejy4bu/backyardEvolution/metadata/testSamples.txt)%10 ~/drosophila_polymorph/backyardEvol/3.markdup.sh
# sbatch --array=1-$(wc -l < /scratch/ejy4bu/backyardEvolution/metadata/allSamples.txt)%20 ~/drosophila_polymorph/backyardEvol/3.markdup.sh
# SAMPLE_LIST="/scratch/ejy4bu/backyardEvolution/metadata/allSamples.txt"
SAMPLE_LIST="/scratch/ejy4bu/backyardEvolution/metadata/testSamples.txt"

sampName=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
SAMPLE_DIR="/scratch/ejy4bu/backyardEvolution/fastq/${sampName}"

# # to test: Dsimu_m_albe_2020_11_01_0092
# SAMPLE_DIR=/scratch/ejy4bu/backyardEvolution/fastq/Dsimu_m_albe_2020_11_01_0092/
# sampName=$(basename "$SAMPLE_DIR")

module load picard


if [ ! -f "${SAMPLE_DIR}/${sampName}.sorted.bam" ]; then
    echo "no sorted bam. "
    exit 1
fi

if [ ! -f "${SAMPLE_DIR}/${sampName}.markdup.bam" ]; then
    echo "dedup sample ${sampName}. "
    java -Xmx45G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        I="${SAMPLE_DIR}/${sampName}.sorted.bam"  \
        O="${SAMPLE_DIR}/${sampName}.markdup.bam" \
        M="${SAMPLE_DIR}/${sampName}.markdup.metrics.txt" \
        CREATE_INDEX=true

    samtools flagstat "${SAMPLE_DIR}/${sampName}.markdup.bam"
else 
    echo "already markdup'd " 
    exit 0
fi

echo "complete"