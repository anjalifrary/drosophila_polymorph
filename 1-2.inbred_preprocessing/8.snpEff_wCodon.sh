#!/usr/bin/env bash

#SBATCH -J snpeff               # A single job name for the array
#SBATCH --ntasks-per-node=10    # one core
#SBATCH -N 1                    # on one node
#SBATCH -t 0-10:00              # 10 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/snpeff.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/snpeff.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

set -euo pipefail

module load bcftools
module load htslib
module load gcc/11.4.0

SNPEFF=/project/berglandlab/multispecies_endemism/snpEFF/v4.3t/snpEff/

### MEL DGRP2:
# outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/"
# in_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.vcf.gz"
# out_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.vcf.gz"

# echo "Annotating Dm6 vcf with SnpEff..."

# java -Xmx32G \
#     -jar ${SNPEFF}/snpEff.jar ann \
#     -formatEff \
#     -v BDGP6.86 \
#     -stats ${outdir}/snpEff_summary.html \
#     $in_vcf \
#     | bgzip -@ 10 -c - > \
#     $out_vcf

# ### SIM Signor:
# outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/"
# in_vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.vcf.gz"
# out_vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.vcf.gz"

# echo "Annotating Dsim3 vcf with SnpEff..."

# java -Xmx32G \
#     -jar ${SNPEFF}/snpEff.jar ann \
#     -formatEff \
#     Dsim_v3.1 \
#     -stats ${outdir}/snpEff_summary.html \
#     $in_vcf \
#     | bgzip -@ 10 -c - > \
#     $out_vcf

#########################################

echo "Finished annotating"
echo "Indexing..."

bcftools index $out_vcf

echo "Complete"

cp "${out_vcf}" \
/project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/

cp "${out_vcf}.csi" \
/project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/
