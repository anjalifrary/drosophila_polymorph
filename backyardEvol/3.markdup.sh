#!/usr/bin/env bash
#
#SBATCH -J be_dup # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/ejy4bu/err_outs/be/markdup.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/be/markdup.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


# module load gcc/11.4.0 htslib
# module load sratoolkit/3.1.1
# module load bwa
# module load samtools
module load picard
# module load gatk
# module load fastqc

set -euo pipefail

### to Run
# sbatch --array=1-$( wc -l < /scratch/ejy4bu/backyardEvolution/metadata/sim_samps.csv )%10 ~/drosophila_polymorph/backyardEvol/1.removeAdapters.sh
# sbatch --array=1-$( wc -l < /scratch/ejy4bu/backyardEvolution/metadata/mel_samps.csv )%10 ~/drosophila_polymorph/backyardEvol/1.removeAdapters.sh

# ref_sim="/scratch/ejy4bu/backyardEvolution/references/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"
# ref_mel="/scratch/ejy4bu/backyardEvolution/references/dmel-all-chromosome-r6.12.fasta"

# to test: Dsimu_m_albe_2020_11_01_0092

SAMPLE_DIR=/scratch/ejy4bu/backyardEvolution/fastq/Dsimu_m_albe_2020_11_01_0092/
sampName=$(basename $SAMPLE_DIR)

### 4. dedup mark duplicates (GATK)


if [ ! -f "${SAMPLE_DIR}/${sampName}.sorted.bam" ]; then
    echo "no sorted bam. "
    exit
fi

if [ ! -f "${SAMPLE_DIR}/${sampName}.markdup.bam" ]; then
    echo "dedup sample ${sampName}. "
    java -Xmx45G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        I="${SAMPLE_DIR}/${sampName}.sorted.bam"  \
        O="${SAMPLE_DIR}/${sampName}.markdup.bam" \
        M="${SAMPLE_DIR}/${sampName}.markdup.metrics.txt" \
        CREATE_INDEX=true

    samtools flagstat "${SAMPLE_DIR}/${sampName}.markdup.bam"
fi

echo "complete"