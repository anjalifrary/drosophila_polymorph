#!/bin/bash
#
#SBATCH -J be_checkbam # A single job name for the array
#SBATCH --cpus-per-task=1
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00 # 10 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/ejy4bu/err_outs/be/check/checkBAM.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/be/check/checkBAM.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

set -euo pipefail
module load samtools


SAMPLE_LIST="/scratch/ejy4bu/backyardEvolution/metadata/allSamples.txt"
wd="/scratch/ejy4bu/backyardEvolution/fastq"


missing=0
empty=0
corrupt=0
missingIndex=0
success=0

while IFS= read -r sample; do 
    SAMPLE_DIR="${wd}/${sample}"
    bam="${wd}/${sample}/${sample}.sorted.markdup.clipped.bam"
    if [ ! -f "$bam" ]; then
        echo "MISSING: $bam"
        missing=$((missing+1))
        continue
    fi
    if [ ! -s "$bam" ]; then
        echo "EMPTY:   $bam"
        empty=$((empty+1))
        continue
    fi
    if ! samtools quickcheck "$bam"; then
        echo "CORRUPT: $bam"
        corrupt=$((corrupt+1))
        continue
    fi
    if [ ! -f "${bam}.bai" ]; then
        echo "MISSING_INDEX: $bam"
        missingIndex=$((missingIndex+1))
        continue
    fi
    success=$((success+1))
    echo "PASSED: $sample"

    echo "Removing intermediates for $sample..."

    rm -f "${wd}/${sample}/${sample}.sorted.bam"
    rm -f "${wd}/${sample}/${sample}.sorted.bam.bai"
    rm -f "${wd}/${sample}/${sample}.sorted.markdup.bam"
    rm -f "${wd}/${sample}/${sample}.sorted.markdup.bai"
    rm -f "${wd}/${sample}/${sample}.trimmed_1.fq.gz"
    rm -f "${wd}/${sample}/${sample}.trimmed_2.fq.gz"
    rm -f "${wd}/${sample}/${sample}_1.fq.gz"
    rm -f "${wd}/${sample}/${sample}_2.fq.gz"

    echo "Cleanup complete: $sample"


done < "$SAMPLE_LIST"

echo
echo "=============================="
echo "Successful samples: $success"
echo "Missing files:      $missing"
echo "Empty files:        $empty"
echo "Corrupt files:      $corrupt"
echo "Missing index:      $missingIndex"
echo "=============================="

if (( $missing + $empty + $corrupt + $missingIndex == 0 )); then
    echo "All samples passed: $success ..  Copying BAMs..."
    new_dir="/project/berglandlab/anjali/be_flies/bam/"

    while IFS= read -r sample; do
        samp_dir="${new_dir}/${sample}/"
        mkdir -p "${samp_dir}"
        echo "copying $sample..."
        cp "${wd}/${sample}/${sample}.sorted.markdup.clipped.bam" "${samp_dir}"
        cp "${wd}/${sample}/${sample}.sorted.markdup.clipped.bam.bai" "${samp_dir}"
        cp "${wd}/${sample}/${sample}.fastp.json" "${samp_dir}"
        cp "${wd}/${sample}/${sample}.fastp.html" "${samp_dir}"
        cp "${wd}/${sample}/${sample}.sorted.markdup.metrics.txt" "${samp_dir}"

    done < "$SAMPLE_LIST"

    echo "Copy complete."
else
    echo "Errors detected. BAMs were not copied."
    exit 1
fi

# ### get missing samples
# SAMPLE_LIST="/scratch/ejy4bu/backyardEvolution/metadata/allSamples.txt"
# wd="/scratch/ejy4bu/backyardEvolution/fastq"

# > /scratch/ejy4bu/backyardEvolution/metadata/missingSamples.txt

# while IFS= read -r sample; do
#     bam="${wd}/${sample}/${sample}.sorted.markdup.clipped.bam"
#     if [ ! -f "$bam" ]; then
#         echo "$sample" >> /scratch/ejy4bu/backyardEvolution/metadata/missingSamples.txt
#     fi
# done < "$SAMPLE_LIST"

# ### remove og fastq files to save space
# SAMPLE_LIST="/scratch/ejy4bu/backyardEvolution/metadata/allSamples.txt"
# wd="/scratch/ejy4bu/backyardEvolution/fastq"

# > /scratch/ejy4bu/backyardEvolution/delete_fastqs.txt

# while IFS= read -r sample; do
#     bam="${wd}/${sample}/${sample}.sorted.markdup.clipped.bam"

#     if [ -f "$bam" ] && [ -s "$bam" ] && samtools quickcheck "$bam" && [ -f "${bam}.bai" ]; then
#         echo "${wd}/${sample}/${sample}_1.fq.gz" >> /scratch/ejy4bu/backyardEvolution/delete_fastqs.txt
#         echo "${wd}/${sample}/${sample}_2.fq.gz" >> /scratch/ejy4bu/backyardEvolution/delete_fastqs.txt
#     fi
# done < "$SAMPLE_LIST"

# echo "Files marked for deletion:"
# cat /scratch/ejy4bu/backyardEvolution/delete_fastqs.txt
# echo
# echo "Number of files: $(wc -l < /scratch/ejy4bu/backyardEvolution/delete_fastqs.txt)"

# less /scratch/ejy4bu/backyardEvolution/delete_fastqs.txt

# # to delete:
# xargs -d '\n' rm -v < /scratch/ejy4bu/backyardEvolution/delete_fastqs.txt

# hdquota -s
