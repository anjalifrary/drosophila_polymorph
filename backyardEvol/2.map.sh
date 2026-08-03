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

### to Run
# sbatch --array=1-$( wc -l < /scratch/ejy4bu/backyardEvolution/metadata/sim_samps.csv )%10 ~/drosophila_polymorph/backyardEvol/1.removeAdapters.sh
# sbatch --array=1-$( wc -l < /scratch/ejy4bu/backyardEvolution/metadata/mel_samps.csv )%10 ~/drosophila_polymorph/backyardEvol/1.removeAdapters.sh


# to test: Dsimu_m_albe_2020_11_01_0092



SAMPLE_DIR=/scratch/ejy4bu/backyardEvolution/fastq/Dsimu_m_albe_2020_11_01_0092/
sampName=$(basename $SAMPLE_DIR)


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


bwa mem \
    -t 8 \
    reference.fa \
    sample_R1.trim.fastq.gz \
    sample_R2.trim.fastq.gz |
samtools view -b |
samtools sort -@8 -o sample.sorted.bam

samtools index sample.sorted.bam