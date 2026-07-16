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
module load htslib
module load gcc/11.4.0

SNPEFF=/project/berglandlab/multispecies_endemism/snpEFF/v4.3t/snpEff/

outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/"

# not reheadered:
# cp /project/berglandlab/Dmel_genomic_resources/DGRP/vcf/DGRP2.source_BCM-HGSC.dm6.final.vcf.gz /scratch/ejy4bu/drosophila/inbred/combined_vcf/
in_vcf="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2.source_BCM-HGSC.dm6.final.vcf.gz"


# in_vcf="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2.source_BCM-HGSC.dm6.final.newheader.vcf.gz"
out_vcf="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2.source_BCM-HGSC.dm6.final.newheader.ann.eff.vcf.gz"

# file /scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2.source_BCM-HGSC.dm6.final.newheader.vcf.gz.gz
# bcftools view -h /scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2.source_BCM-HGSC.dm6.final.newheader.vcf.gz.gz | head
# file /project/berglandlab/Dmel_genomic_resources/DGRP/vcf/DGRP2.source_BCM-HGSC.dm6.final.vcf.gz 
# bcftools view -h /project/berglandlab/Dmel_genomic_resources/DGRP/vcf/DGRP2.source_BCM-HGSC.dm6.final.vcf.gz | head
# file $in_vcf
# bcftools view -h $in_vcf | head

# file /project/berglandlab/DGRP_freeze2_vcf/DGRP2.source_BCM-HGSC.dm6.final.newheader.vcf.gz

echo "Annotating Dm6 vcf with SnpEff..."

java -Xmx32G \
    -jar ${SNPEFF}/snpEff.jar ann \
    -formatEff \
    -v BDGP6.86 \
    -stats ${outdir}/snpEff_summary.html \
    $in_vcf \
    | bgzip -@ 10 -c - > \
    $out_vcf

echo "Finished annotating"
echo "Indexing..."

bcftools index /scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2.source_BCM-HGSC.dm6.final.newheader.ann.eff.vcf.gz

echo "Complete"

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