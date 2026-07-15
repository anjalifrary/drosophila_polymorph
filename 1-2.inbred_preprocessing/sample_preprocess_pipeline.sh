#!/usr/bin/env bash
#
#SBATCH -J preprocess # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/prepipeline.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/prepipeline.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=2-183

set -euo pipefail


#ijob -A berglandlab -c10 -p standard --mem=100G

ref="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"


# ### do once
# bwa index ${ref}
# samtools faidx ${ref}
# gatk CreateSequenceDictionary -R ${ref} 

module load gcc/11.4.0 htslib
module load sratoolkit/3.1.1
module load bwa
module load samtools
module load picard
module load gatk
module load fastqc



MY_DATA="/scratch/ejy4bu/drosophila/inbred/fastq/PRJNA318623/"
# SAMPLES=($(ls -d ${MY_DATA}*/))
SAMPLES=("${MY_DATA}"*/)

echo "Samples = ${#SAMPLES[@]}"

SAMPLE_DIR="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
samp_name=$(basename "$SAMPLE_DIR")

fastq="${SAMPLE_DIR}/${samp_name}.fastq.gz"


echo "Processing sample ${samp_name}. (Array task ID: $SLURM_ARRAY_TASK_ID)"

mkdir -p ${SAMPLE_DIR}/{fastqc,fastp,bam,gvcf}

### 1. fastqc
fastqc \
    -o ${SAMPLE_DIR}/fastqc \
    -t 10 \
    ${fastq}

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

### 3. bwa mem single end mode -> to sorted bam
if [ ! -f "${SAMPLE_DIR}/bam/${samp_name}.sorted.bam" ]; then
    echo "mapping/sorting sample ${samp_name}. "
    bwa mem \
        -t 10 \
        -K 100000000 \
        -Y \
        -R "@RG\tID:${samp_name}\tSM:${samp_name}\tPL:ILLUMINA\tLB:${samp_name}" \
        ${ref} \
        ${SAMPLE_DIR}/fastp/${samp_name}.trimmed.fastq.gz \
        | samtools view -uh -q 20 \
        | samtools sort --threads 2 -o ${SAMPLE_DIR}/bam/${samp_name}.sorted.bam -

    samtools index ${SAMPLE_DIR}/bam/${samp_name}.sorted.bam
    echo "finished mapping $samp_name"
fi

### 4. dedup mark duplicates (GATK)

if [ ! -f "${SAMPLE_DIR}/bam/${samp_name}.markdup.bam" ]; then
    echo "dedup sample ${samp_name}. "
    java -Xmx45G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        I=${SAMPLE_DIR}/bam/${samp_name}.sorted.bam  \
        O=${SAMPLE_DIR}/bam/${samp_name}.markdup.bam \
        M=${SAMPLE_DIR}/bam/${samp_name}.markdup.metrics.txt \
        CREATE_INDEX=true
fi
# should i remove duplicates or just flag?

### 5. Haplotype caller 
if [ ! -f "${SAMPLE_DIR}/gvcf/${samp_name}.g.vcf.gz" ]; then
    echo "haplotype caller on sample ${samp_name}. "
    gatk --java-options "-Xmx40G -XX:+UseG1GC" HaplotypeCaller \
    -R ${ref} \
    -I ${SAMPLE_DIR}/bam/${samp_name}.markdup.bam \
    -O ${SAMPLE_DIR}/gvcf/${samp_name}.g.vcf.gz \
    -ERC GVCF \
    --native-pair-hmm-threads 8
fi

echo "Complete"

# ### 6. concatenate vcfs 
# module load bcftools

# # Define input and output directories
# VCF_DIR="/project/berglandlab/chlorella_sequencing/vcfs"
# OUTPUT_VCF="/project/berglandlab/chlorella_sequencing/combined_chlorella.vcf.gz"

# # Create a list of VCF files
# VCF_LIST=$(ls ${VCF_DIR}/*.vcf.gz | sort)

# # Concatenate VCFs using bcftools
# bcftools concat -a -O z -o ${OUTPUT_VCF} ${VCF_DIR}/*.vcf.gz

# # Index the final VCF (optional but recommended for SNPeff)
# bcftools index ${OUTPUT_VCF}

# ### 7. SnpEFF
# SNPEFF_DIR=$HOME/snpEff
# VCF_IN="/project/berglandlab/chlorella_sequencing/combined_chlorella.vcf.gz"
# VCF_OUT="/project/berglandlab/chlorella_sequencing/combined_chlorella_annotated.vcf"

# java -Xmx32g -jar $SNPEFF_DIR/snpEff.jar -v chlorella $VCF_IN > $VCF_OUT