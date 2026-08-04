#!/usr/bin/env bash
#
#SBATCH -J be_pipeEach # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-18:00 # 18 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/ejy4bu/err_outs/be/pipe/pipeEach.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/be/pipe/pipeEach.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=501-1000%20

# 1-500
# 501-1000
# 1001-1500
# 1501-2000
# 2001-2500
# 2501-3000
# 3001-3500
# 3501-3665


set -euo pipefail

SAMPLE_LIST="/scratch/ejy4bu/backyardEvolution/metadata/allSamples.txt"
# SAMPLE_LIST="/scratch/ejy4bu/backyardEvolution/metadata/testSamples.txt"

sampName=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
SAMPLE_DIR="/scratch/ejy4bu/backyardEvolution/fastq/${sampName}"
echo "Processing ${sampName}"


# 1. remove adapters (trim)

if [ ! -f "${SAMPLE_DIR}/${sampName}_1.fq.gz" ]; then
    echo "no fastq 1 file. "
    exit 1
fi

if [ ! -f "${SAMPLE_DIR}/${sampName}_2.fq.gz" ]; then
    echo "no fastq 2 file. "
    exit 1
fi

if [ ! -f "${SAMPLE_DIR}/${sampName}.trimmed_1.fq.gz" ]; then 
    module load miniforge && conda activate fastp
    echo "trimming sample ${sampName}. "
    fastp \
        --in1 "${SAMPLE_DIR}/${sampName}_1.fq.gz" \
        --in2 "${SAMPLE_DIR}/${sampName}_2.fq.gz" \
        --out1 "${SAMPLE_DIR}/${sampName}.trimmed_1.fq.gz" \
        --out2 "${SAMPLE_DIR}/${sampName}.trimmed_2.fq.gz" \
        --json "${SAMPLE_DIR}/${sampName}.fastp.json" \
        --html "${SAMPLE_DIR}/${sampName}.fastp.html" \
        --detect_adapter_for_pe \
        --thread ${SLURM_CPUS_PER_TASK}

    conda deactivate
    conda deactivate
else 
    echo "already trimmed " 
fi

echo "trimming complete..."

# 2. mapping

module load bwa
module load samtools

# indexed ref must exist in same directory
ref_sim="/scratch/ejy4bu/backyardEvolution/references/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"
ref_mel="/scratch/ejy4bu/backyardEvolution/references/dmel-all-chromosome-r6.12.fasta"

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
    samtools sort --threads ${SLURM_CPUS_PER_TASK} -o ${SAMPLE_DIR}/${sampName}.sorted.bam -

    echo "mapped and sorted"
    samtools index ${SAMPLE_DIR}/${sampName}.sorted.bam
    echo "indexed"
    samtools flagstat ${SAMPLE_DIR}/${sampName}.sorted.bam

    if ! samtools quickcheck "${SAMPLE_DIR}/${sampName}.sorted.bam"; then
        echo "ERROR: corrupt sorted bam"
        exit 1
    fi
else 
    echo "already mapped & sorted " 
fi

echo "mapping/index complete"

# 3. mark duplicates

module load picard
module load samtools


if [ ! -f "${SAMPLE_DIR}/${sampName}.sorted.bam" ]; then
    echo "no sorted bam. "
    exit 1
fi

if [ ! -f "${SAMPLE_DIR}/${sampName}.markdup.bam" ]; then
    echo "mark dup sample ${sampName}. "
    java -Xmx45G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        I="${SAMPLE_DIR}/${sampName}.sorted.bam"  \
        O="${SAMPLE_DIR}/${sampName}.markdup.bam" \
        M="${SAMPLE_DIR}/${sampName}.markdup.metrics.txt" \
        CREATE_INDEX=true

    samtools flagstat "${SAMPLE_DIR}/${sampName}.markdup.bam"
    
    if ! samtools quickcheck "${SAMPLE_DIR}/${sampName}.markdup.bam"; then
        echo "ERROR: corrupt markdup bam"
        exit 1
    fi
else 
    echo "already markdup'd "
fi

echo "mark dupes complete"

# 4. soft clipping overlapping reads  

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

    if ! samtools quickcheck "${SAMPLE_DIR}/${sampName}.clipped.bam"; then
        echo "ERROR: corrupt clipped bam"
        exit 1
    fi
else 
    echo "already clipped " 
fi

echo "complete"

# all done

echo "all steps complete. bam located at ${SAMPLE_DIR}/${sampName}.clipped.bam"