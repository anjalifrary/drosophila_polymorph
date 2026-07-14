#!/usr/bin/env bash
#
#SBATCH -J nextflow_sarek # job name
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/sarek.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/sarek.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module purge
module load java/17
module load nextflow
module load apptainer

ref="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"
sample_csv="/scratch/ejy4bu/drosophila/inbred/sarek/sample_inputs_sarek.csv"
out_dir="/scratch/ejy4bu/drosophila/inbred/sarek_results/"

export NXF_APPTAINER_CACHEDIR="/scratch/ejy4bu/tmp/apptainer_cache"
mkdir -p $NXF_APPTAINER_CACHEDIR

nextflow run nf-core/sarek \
    -r 3.4.0 \
    -profile apptainer \
    --input $sample_csv \
    --outdir $out_dir \
    --fasta $ref \
    --save_reference \
    --aligner bwa-mem \
    --tools haplotypecaller \
    --igenomes_ignsore \
    -resume

# set --tools haplotypecaller restricts it to GATK Haplotype Caller 
# aligner as bwa-mem skips other defaults (Dragmap or bwa-mem2)
# igenomes_ignore to ignore built in iGenomes paths for humans / mouse