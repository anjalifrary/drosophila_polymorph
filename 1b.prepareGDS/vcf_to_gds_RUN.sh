#!/usr/bin/env bash
#SBATCH -J vcf2gds         # Job name
#SBATCH --ntasks=1         # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1               # Run on one node
#SBATCH -t 0-10:00         # 10 hours runtime
#SBATCH --mem=100G         # Memory per node
#SBATCH -o /scratch/ejy4bu/err_outs/GDS/vcf2gds.%A_%a.out  # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/GDS/vcf2gds.%A_%a.err  # Standard error
#SBATCH -p standard        # Partition
#SBATCH --account=berglandlab

export R_LIBS_USER=~/Rlibs
module load gcc/11.4.0 openmpi/4.1.4 
module load R/4.5.0
module load bcftools

# reheader dmel file
vcf="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/DGRP2.source_BCM-HGSC.dm6.final.norm.filtered.ann.eff.vcf.gz"
# bcftools view -h $vcf \
# > /scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/DRGP2.dm6.header.txt

# old: ##INFO=<ID=NS,Number=0,Type=Flag,Description="Non-synonymous SNP">
# new: ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">

bcftools reheader \
    -h /scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/DRGP2.dm6.header.txt \
    -o /scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/DGRP2.source_BCM-HGSC.dm6.final.norm.reheadered.ann.eff.vcf.gz \
    $vcf

bcftools index /scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/DGRP2.source_BCM-HGSC.dm6.final.norm.reheadered.ann.eff.vcf.gz

Rscript 1b.prepareGDS/vcf_to_gds.R

# cp /scratch/ejy4bu/drosophila/gds_files/DGRP2.source_BCM-HGSC.dm6.final.norm.ann.eff.gds \
# /project/berglandlab/anjali/drosophila_polymorphism/data_files/gds/