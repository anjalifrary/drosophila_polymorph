#!/usr/bin/env bash
#
#SBATCH -J gdbi # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/gdbi.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/gdbi.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
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

#  /project/berglandlab/alan/privatePolymorphisms/simulans/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.vcf.gz | \

chroms=("2L" "2R" "3L" "3R" "4" "X")
chr=${chroms[$SLURM_ARRAY_TASK_ID]}

JAVAMEM=80G
CPU=10

echo "Processing Chromosome: ${chr}"

# module load sratoolkit/3.1.1
# module load bwa
# module load samtools
# module load picard
# module load gatk
# module load fastqc

### 1. get gVCFs
# # find all vcfs
# find ${gvcf_dir} -name "*.g.vcf.gz" \
#     | sort > ${outdir}/gvcf.list

# echo "Found $(wc -l < ${outdir}/gvcf.list) gVCFs"
# if [ "$(wc -l < ${outdir}/gvcf.list)" -eq 0 ]; then
#     echo "ERROR: No gVCFs found in ${gvcf_dir}!" >&2
#     exit 1
# fi

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

### from signor

# /home/cmb-07/sn1/ssignor/jre1.7.0/bin/java -Xmx20g -jar /home/cmb-07/sn1/ssignor/next-gen_example_sim/bin/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
# 	-T GenotypeGVCFs \
# 	-R $REF_PFX \
# 	-L X \
#   [the list of --variant ..]
#   -o .vcf

# gatk CombineGVCFs \
#     -R ${ref} \
#     $(while read gvcf; do echo "--variant $gvcf"; done < ${outdir}/gvcf.list) \
#     -O ${outdir}/dsim3.signor.combined.g.vcf.gz
# echo "Done combining gVCFs"

### Consolidate vcfs 

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
    --batch-size 50 \
    --sample-name-map ${sample_map} \
    --reader-threads $CPU \
    -L ${chr}

rm -rf "$tmp"

echo "genomicsdbimport ${chr} completed successfully!"