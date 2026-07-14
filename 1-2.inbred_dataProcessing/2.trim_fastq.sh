#!/usr/bin/env bash
#
#SBATCH -J trimSRA # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/trim.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/trim.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-3
#ijob -A berglandlab -c10 -p standard --mem=40G


# Load necessary modules
module load gcc htslib
module load trimmomatic

#for pipeline, get sample_dir from array
wd="/scratch/ejy4bu/drosophila/inbred/fastq/PRJNA318623"
# sample_dir="$1"
# sample_dir="${wd}/SRR3585384"
sample_dir=$(find "$wd" -maxdepth 1 -type d -name "SRR*" | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "processing ${sample_dir}"

samp_name=$(basename "$sample_dir")

if [ -f "${sample_dir}/${samp_name}.fastq" ]; then
  gzip ${sample_dir}/${samp_name}.fastq
fi

input="${sample_dir}/${samp_name}.fastq.gz"
output="${sample_dir}/${samp_name}.trimmed.fastq.gz"

if [ -f "$output" ]; then
    echo "Already trimmed: ${samp_name}"
    exit 0
fi

echo "Running Trimmomatic on ${samp_name}"
# note: running trimmomatic on SE (single-end) mode & not running PEAR anymore (single-end reads)

java -jar $EBROOTTRIMMOMATIC/trimmomatic-*.jar SE \
    -threads 10 \
    "$input" \
    "$output" \
    ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-SE.fa:2:30:10

echo "Done trimming ${samp_name}"


### questions:
#  trimming parameters? quality params?
    # LEADING:3 \
    # TRAILING:3 \
    # SLIDINGWINDOW:4:20 \
    # MINLEN:36


for samp in ${sample_dir}/${samp_name}.fastq.gz; do
    if ls "$sample_dir"/*.trimm.fastq 1> /dev/null 2>&1; then
        echo "sample already trimmed: $samp_name"
        continue
    fi

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-*.jar SE -threads 10 \
        ${lane_name}_1.fq.gz \
        ${lane_name}_2.fq.gz \
        ${lane_name}_1.P.trimm.fastq \
        ${lane_name}_1.U.trimm.fastq \
        ${lane_name}_2.P.trimm.fastq \
        ${lane_name}_2.U.trimm.fastq \
        ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:8:true

    /home/ejy4bu/miniconda3/bin/pear \
        -f ${lane_name}_1.P.trimm.fastq \
        -r ${lane_name}_2.P.trimm.fastq \
        -o ${lane_name} \
        -j 10
    echo "Done working on $lane_name"

done
echo "All lanes trimmed for $samp_name"






for forward in "$sample_dir"/*_1.fq.gz; do
    reverse="${forward/_1.fq.gz/_2.fq.gz}"
    lane_name="${forward%_1.fq.gz}"

    if ls "$sample_dir"/*_1.P.trimm.fastq 1> /dev/null 2>&1; then
        echo "sample already trimmed: $samp_name"
        continue
    fi

    echo "Processing lane: ${lane_name}"

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-*.jar PE -threads 10 \
        ${lane_name}_1.fq.gz \
        ${lane_name}_2.fq.gz \
        ${lane_name}_1.P.trimm.fastq \
        ${lane_name}_1.U.trimm.fastq \
        ${lane_name}_2.P.trimm.fastq \
        ${lane_name}_2.U.trimm.fastq \
        ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:8:true

    /home/ejy4bu/miniconda3/bin/pear \
        -f ${lane_name}_1.P.trimm.fastq \
        -r ${lane_name}_2.P.trimm.fastq \
        -o ${lane_name} \
        -j 10
    echo "Done working on $lane_name"

done
echo "All lanes trimmed for $samp_name"

