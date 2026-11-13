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

#########################################
## D. melanogaster:
outdir_mel="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/filtering/mel/"
in_mel="${outdir_mel}/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.snpGap10.snpsOnly.repeatmasked.wmdust.vcf.gz"
out_mel="${outdir_mel}/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.snpGap10.snpsOnly.repeatmasked.wmdust.ann.eff.vcf.gz"

echo "Annotating Dm6 vcf with SnpEff..."

java -Xmx32G \
    -jar ${SNPEFF}/snpEff.jar ann \
    -formatEff \
    -v BDGP6.86 \
    -stats ${outdir_mel}/snpEff_summary.mel.html \
    $in_mel \
    | bgzip -@ 10 -c - > \
    $out_mel

echo "Finished annotating"
echo "Indexing..."

bcftools index $out_mel

echo "D. mel complete"

#########################################
### D. simulans:
outdir_sim="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/filtering/sim/"
in_sim="${outdir_sim}/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.norm.snpGap10.snpsOnly.repeatmasked.wmdust.vcf.gz"
out_sim="${outdir_sim}/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.norm.snpGap10.snpsOnly.repeatmasked.wmdust.ann.eff.vcf.gz"

echo "Annotating Dsim3 vcf with SnpEff..."

java -Xmx32G \
    -jar ${SNPEFF}/snpEff.jar ann \
    -formatEff \
    Dsim_v3.1 \
    -stats ${outdir_sim}/snpEff_summary.html \
    $in_sim \
    | bgzip -@ 10 -c - > \
    $out_sim


echo "Finished annotating"
echo "Indexing..."

bcftools index $out_sim

echo "D. sim complete"
