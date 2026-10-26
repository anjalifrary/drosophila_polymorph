#!/usr/bin/env bash
#
#SBATCH -J gdbi # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-18:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/gdbi.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/gdbi.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --array=0-5

set -euo pipefail

module load gcc/11.4.0
module load gatk
module load bcftools
module load htslib

ref="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"
outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/"
mkdir -p ${outdir}
gvcf_dir="/scratch/ejy4bu/drosophila/inbred/fastq/PRJNA318623"

chroms=("sim_2L" "sim_2R" "sim_3L" "sim_3R" "sim_4" "sim_X")
# chroms=("2L" "2R" "3L" "3R" "4" "X")
chr=${chroms[$SLURM_ARRAY_TASK_ID]}

JAVAMEM=80G
CPU=10

echo "Processing Chromosome: ${chr}"


###### sample map
sample_map="${outdir}/dsim_sample_map.txt"
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    echo "Finding gVCFs and creating sample map..."
    find ${gvcf_dir} -name "*.g.vcf.gz" | sort > ${outdir}/gvcf.list
    
    # Generate GATK-compliant sample map: [SampleName]\t[Path]
    awk '{n=$0; sub(".*/", "", n); sub(/\..*/, "", n); print n "\t" $0}' \
        ${outdir}/gvcf.list > ${sample_map}
        
    echo "Found $(wc -l < ${sample_map}) gVCFs."
    else
    echo "Task ${SLURM_ARRAY_TASK_ID} waiting for Task 0 to generate sample map..."
    while [ ! -f "${sample_map}" ] || [ ! -s "${sample_map}" ]; do
        sleep 5
    done
fi

work="${outdir}/dsim_genomicsdb_${chr}"
tmp="/scratch/ejy4bu/tmp/temp_dsim_gdbi_${chr}"
# error handling:
if [ -d "$work" ]; then
    echo "Removing existing GenomicsDB workspace..."
    rm -rf $work
fi

mkdir -p $tmp

echo "Running GenomicsDBImport..."

gatk --java-options "-Xmx${JAVAMEM}" GenomicsDBImport \
    -R ${ref} \
    --genomicsdb-workspace-path $work \
    --tmp-dir $tmp \
    --batch-size 150 \
    --sample-name-map ${sample_map} \
    --reader-threads $CPU \
    -L ${chr}

rm -rf "$tmp"

echo "genomicsdbimport ${chr} completed successfully."