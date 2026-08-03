#!/usr/bin/env bash
#
#SBATCH -J be_trim # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/ejy4bu/err_outs/be/trim.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/be/trim.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### to Run
# sbatch --array=1-$( wc -l < /scratch/ejy4bu/backyardEvolution/metadata/sim_samps.csv )%10 ~/drosophila_polymorph/backyardEvol/1.removeAdapters.sh
# sbatch --array=1-$( wc -l < /scratch/ejy4bu/backyardEvolution/metadata/mel_samps.csv )%10 ~/drosophila_polymorph/backyardEvol/1.removeAdapters.sh

# to test: Dsimu_m_albe_2020_11_01_0092


module load fastp

SAMPLE_DIR=/scratch/ejy4bu/backyardEvolution/fastq/Dsimu_m_albe_2020_11_01_0092/
sampName=$(basename $SAMPLE_DIR)

if [ ! -f "${SAMPLE_DIR}/${sampName}_1.fq.gz" ]; then
    echo "no fastq 1 file. "
    exit
fi

if [ ! -f "${SAMPLE_DIR}/${sampName}_2.fq.gz" ]; then
    echo "no fastq 2 file. "
    exit
fi

if [ ! -f "${SAMPLE_DIR}/${sampName}.trimmed_1.fq.gz" ]; then 
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
fi

# https://github.com/OpenGene/fastp#adapters
### flag detect_adapter_for_pe 

## For PE data, the auto adapter detection is disabled by default. You can enable it by specifing 
# -2 or --detect_adapter_for_pe. If you want to obtain ultra-clean data, this option is strongly suggested.


## For PE data, the adapters can be trimmed automatically by per-read overlap analysis, 
# which seeks for the overlap of each pair of reads. This method is robust and fast, 
# so normally you don't have to input the adapter sequence. But you can still specify the 
# adapter sequences for read1 by --adapter_sequence, and for read2 by --adapter_sequence_r2. 
# In case fastp fails to find an overlap for some pairs (i.e. due to low quality bases), it will 
# use these sequences to trim adapters for read1 and read2 respectively.


## For PE data, fastp will run a little slower if you specify the sequence adapters or enable the adapter auto-detection. 
# But it may result in a slightly cleaner output (usually finds 0.1% to 0.5% more adapters), 
# since the overlap analysis may fail due to sequencing errors.

## For PE data, you can specify --allow_gap_overlap_trimming to allow up to one gap when trim 
# adapters by overlap analysis for PE data. By default no gap is allowed. This may take more 
# time and usually have very limited effect (finds ~0.01% more adapters).

