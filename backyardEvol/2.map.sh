#!/usr/bin/env bash
#
#SBATCH -J be_map # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/ejy4bu/err_outs/be/map.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/be/map.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load bwa
module load samtools

set -euo pipefail

### to Run
# sbatch --array=1-$(wc -l < /scratch/ejy4bu/backyardEvolution/metadata/allSamples.txt)%20 ~/drosophila_polymorph/backyardEvol/1.removeAdapters.sh
SAMPLE_LIST="/scratch/ejy4bu/backyardEvolution/metadata/allSamples.txt"
sampName=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
SAMPLE_DIR="/scratch/ejy4bu/backyardEvolution/fastq/${sampName}"

# # to test: Dsimu_m_albe_2020_11_01_0092
# SAMPLE_DIR=/scratch/ejy4bu/backyardEvolution/fastq/Dsimu_m_albe_2020_11_01_0092/
# sampName=$(basename "$SAMPLE_DIR")

ref_sim="/scratch/ejy4bu/backyardEvolution/references/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"
ref_mel="/scratch/ejy4bu/backyardEvolution/references/dmel-all-chromosome-r6.12.fasta"

# # do once:
# bwa index ${ref_sim}
# bwa index ${ref_mel}

if [ ! -f "${SAMPLE_DIR}/${sampName}.trimmed_1.fq.gz" ]; then
    echo "no trimmed fastq 1 file. "
    exit 1
fi

if [ ! -f "${SAMPLE_DIR}/${sampName}.trimmed_2.fq.gz" ]; then
    echo "no trimmed fastq 2 file. "
    exit 1
fi

if [[ ${sampName} == Dsimu* ]]; then 
    ref=${ref_sim}
elif [[ ${sampName} == Dmela* ]]; then
    ref=${ref_mel}
else 
    echo "unrecognized species " 
    exit 1
fi

if [ ! -f "${SAMPLE_DIR}/${sampName}.sorted.bam" ]; then 
    echo "mapping sample ${sampName}" 
    bwa mem \
        -t ${SLURM_CPUS_PER_TASK} \
        -R "@RG\tID:${sampName}\tSM:${sampName}\tPL:ILLUMINA" \
        ${ref} \
        ${SAMPLE_DIR}/${sampName}.trimmed_1.fq.gz \
        ${SAMPLE_DIR}/${sampName}.trimmed_2.fq.gz |
    samtools sort --threads 2 -o ${SAMPLE_DIR}/${sampName}.sorted.bam -

    echo "mapped and sorted"
    samtools index ${SAMPLE_DIR}/${sampName}.sorted.bam
    echo "indexed"
    samtools flagstat ${SAMPLE_DIR}/${sampName}.sorted.bam
fi

echo "complete"

# ### 3. bwa mem single end mode -> to sorted bam
# if [ ! -f "${SAMPLE_DIR}/bam/${samp_name}.sorted.bam" ]; then
#     echo "mapping/sorting sample ${samp_name}. "
#     bwa mem \
#         -t 10 \
#         -K 100000000 \
#         -Y \
#         -R "@RG\tID:${samp_name}\tSM:${samp_name}\tPL:ILLUMINA\tLB:${samp_name}" \
#         ${ref} \
#         ${SAMPLE_DIR}/fastp/${samp_name}.trimmed.fastq.gz \
#         | samtools view -uh -q 20 \
#         | samtools sort --threads 2 -o ${SAMPLE_DIR}/bam/${samp_name}.sorted.bam -

#     samtools index ${SAMPLE_DIR}/bam/${samp_name}.sorted.bam
#     echo "finished mapping $samp_name"
# fi


# bwa mem \
#     -t 8 \
#     reference.fa \
#     sample_R1.trim.fastq.gz \
#     sample_R2.trim.fastq.gz |
# samtools view -b |
# samtools sort -@8 -o sample.sorted.bam

# samtools index sample.sorted.bam