#!/usr/bin/env bash

#SBATCH -J runFASTQC # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 4:00:00 ### 15 seconds
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/logs/demo_1.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/logs/demo_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run as: sbatch ~/misc/1000G/done/Dsim/5.DSim.snpEff_codon.sh
### sacct -j 10700313
### cat /scratch/aob2x/logs/demo_1.10700296*.err


 SNPEFF_HOME=/project/berglandlab/multispecies_endemism/snpEFF/v4.3t/snpEff

### annotate
#  module load htslib
#
#  java -jar $SNPEFF_HOME/snpEff.jar ann -formatEff \
#  Dsim_v3.1 \
#  /project/berglandlab/alan/privatePolymorphisms/simulans/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.vcf.gz | \
#  bgzip -@ 20 -c - > \
#  /project/berglandlab/alan/privatePolymorphisms/simulans/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.vcf.gz

### convert
#  module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1
#
#  Rscript --vanilla ~/DESTv3/snpCalling_dev/scatter_gather_annotate/vcf2gds.R \
#  /project/berglandlab/alan/privatePolymorphisms/simulans/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.vcf.gz \
#  20

module load htslib

zcat /project/berglandlab/alan/privatePolymorphisms/simulans/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.vcf.gz |
grep -E "#|missense|synonymous" | \
bgzip -@ 10 -c - > /project/berglandlab/multispecies_endemism/data/Dsim/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.NS_SYN.ann.eff.vcf.gz

module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

Rscript --vanilla ~/DESTv3/snpCalling_dev/scatter_gather_annotate/vcf2gds.R \
/project/berglandlab/multispecies_endemism/data/Dsim/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.NS_SYN.ann.eff.vcf.gz \
20