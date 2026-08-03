#!/usr/bin/env bash
#
#SBATCH -J be_clip # A single job name for the array
#SBATCH --cpus-per-task=2
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/ejy4bu/err_outs/be/clip.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/be/clip.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


set -euo pipefail

### to Run
# sbatch --array=1-$(wc -l < /scratch/ejy4bu/backyardEvolution/metadata/allSamples.txt)%20 ~/drosophila_polymorph/backyardEvol/1.removeAdapters.sh
# SAMPLE_LIST="/scratch/ejy4bu/backyardEvolution/metadata/allSamples.txt"
# sampName=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
# SAMPLE_DIR="/scratch/ejy4bu/backyardEvolution/fastq/${sampName}"

# to test: Dsimu_m_albe_2020_11_01_0092
SAMPLE_DIR=/scratch/ejy4bu/backyardEvolution/fastq/Dsimu_m_albe_2020_11_01_0092/
sampName=$(basename "$SAMPLE_DIR")


if [ ! -f "${SAMPLE_DIR}/${sampName}.markdup.bam" ]; then
    echo "no markdup'd bam. "
    exit 1
fi

if [ ! -f "${SAMPLE_DIR}/${sampName}.clipped.bam" ]; then
    module load miniforge && conda activate bamutil
    echo "soft-clip overlapping reads... sample ${sampName}. "
    bam clipOverlap \
        --in "${SAMPLE_DIR}/${sampName}.markdup.bam" \
        --out "${SAMPLE_DIR}/${sampName}.clipped.bam" \
        --stats

    echo "clipping complete"
    conda deactivate
    conda deactivate

    module load samtools

    samtools index "${SAMPLE_DIR}/${sampName}.clipped.bam"
    echo "indexing complete"

    samtools flagstat "${SAMPLE_DIR}/${sampName}.clipped.bam"
fi

echo "complete"

# https://genome.sph.umich.edu/w/index.php?title=BamUtil:_clipOverlap&mobileaction=toggle_view_desktop
    