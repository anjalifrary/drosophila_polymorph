#!/bin/bash
#
#SBATCH -J be_checkbam # A single job name for the array
#SBATCH --cpus-per-task=8
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/ejy4bu/err_outs/be/checkBAM.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/be/checkBAM.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

set -euo pipefail
module load samtools


SAMPLE_LIST="/scratch/ejy4bu/backyardEvolution/metadata/testSamples.txt"
wd="/scratch/ejy4bu/backyardEvolution/fastq"


missing=0
empty=0
success=0

while IFS= read -r sample; do 
    SAMPLE_DIR="${wd}/${sample}"
    # echo "checking $sample"
    bam="${wd}/${sample}/${sample}.clipped.bam"
    if [ ! -f "$bam" ]; then
        echo "MISSING: $bam"
        continue
    fi
    if [ ! -s "$bam" ]; then
        echo "EMPTY:   $bam"
        continue
    fi
    if ! samtools quickcheck "$bam"; then
        echo "CORRUPT: $bam"
        continue
    fi
    if [ ! -f "${bam}.bai" ]; then
        echo "MISSING_INDEX: $bam"
        continue
    fi
    success=$((success+1))

done < "$SAMPLE_LIST"

echo
echo "=============================="
echo "Successful samples: $success"
echo "Missing files:      $missing"
echo "Empty files:        $empty"
echo "=============================="