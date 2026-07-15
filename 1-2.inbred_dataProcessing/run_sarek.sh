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

cd 1-2.inbred_dataProcessing

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
    -r 3.9.0 \
    -profile apptainer \
    --input $sample_csv \
    --outdir $out_dir \
    --genome null \
    --igenomes_ignore \
    --fasta $ref \
    --save_reference \
    --aligner bwa-mem \
    --trim_fastq \
    --step mapping \
    --tools haplotypecaller \
    --save_trimmed \
    --save_mapped \
    --save_output_as_bam 
    
    # \
    # -resume

# set --tools haplotypecaller restricts it to GATK Haplotype Caller 
# aligner as bwa-mem skips other defaults (Dragmap or bwa-mem2)
# igenomes_ignore to ignore built in iGenomes paths for humans / mouse * and genome null 
# default start in pipieline is --step mapping 
# --save-trimmed saves intermediate trimmed files

# other flags
# --length-required for minimum length of reads to keep
#

#FastQC gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the FastQC help pages.
# The FastQC plots displayed in the MultiQC report shows untrimmed reads. They may contain adapter sequence and potentially regions with low quality.

# FastP is a tool designed to provide all-in-one preprocessing for FastQ files and is used for trimming and splitting. The tool then determines QC metrics for the processed reads.

