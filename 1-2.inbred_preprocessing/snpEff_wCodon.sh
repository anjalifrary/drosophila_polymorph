#!/usr/bin/env bash

#SBATCH -J snpeff               # A single job name for the array
#SBATCH --ntasks-per-node=10    # one core
#SBATCH -N 1                    # on one node
#SBATCH -t 0-10:00              # 10 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/snpeff.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/snpeff.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load bcftools

 SNPEFF=/project/berglandlab/multispecies_endemism/snpEFF/v4.3t/snpEff/

# Dmel6 reference genome


 
echo "Annotating Dm6 vcf with SnpEff..."

java -Xmx32G \
    -jar ${SNPEFF}/snpEff.jar ann -formatEff \
    -v BDGP6.86 \
    -stats ${outdir}/snpEff_summary.html \
    /scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2.source_BCM-HGSC.dm6.final.newheader.vcf.gz \
    | bgzip -@ 10 -c - > \
    /scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2.source_BCM-HGSC.dm6.final.newheader.ann.eff.vcf.gz

bcftools index /scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2.source_BCM-HGSC.dm6.final.newheader.ann.eff.vcf.gz




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

# module load htslib

# zcat /project/berglandlab/alan/privatePolymorphisms/simulans/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.vcf.gz |
# grep -E "#|missense|synonymous" | \
# bgzip -@ 10 -c - > /project/berglandlab/multispecies_endemism/data/Dsim/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.NS_SYN.ann.eff.vcf.gz

# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

# Rscript --vanilla ~/DESTv3/snpCalling_dev/scatter_gather_annotate/vcf2gds.R \
# /project/berglandlab/multispecies_endemism/data/Dsim/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.NS_SYN.ann.eff.vcf.gz \
# 20