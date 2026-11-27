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

species="mel"       # mel or sim

if [[ "$species" == "mel" ]]; then
    outdir="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/filtering/mel/"
    infile="${outdir}/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.snpGap10.snpsOnly.repeatmasked.wmdust.vcf.gz"
    outfile="${outdir}/testGeneID/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.norm.snpGap10.snpsOnly.repeatmasked.wmdust.ann.eff.vcf.gz"

    echo "Annotating Dm6 vcf with SnpEff..."

    java -Xmx32G \
        -jar ${SNPEFF}/snpEff.jar ann \
        -formatEff \
        -v BDGP6.86 \
        -geneId \
        -stats ${outdir}/snpEff_summary.mel.html \
        $infile \
        | bgzip -@ 10 -c - > \
        $outfile

    echo "Mel file annotated"

elif [[ "$species" == "sim" ]]; then
    outdir="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/filtering/sim/"
    infile="${outdir}/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.norm.snpGap10.snpsOnly.repeatmasked.wmdust.vcf.gz"
    outfile="${outdir}/testGeneID/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.norm.snpGap10.snpsOnly.repeatmasked.wmdust.ann.eff.vcf.gz"

    echo "Annotating Dsim3 vcf with SnpEff..."

    java -Xmx32G \
        -jar ${SNPEFF}/snpEff.jar ann \
        -formatEff \
        -v Dsim_v3.1 \
        -geneId \
        -stats ${outdir}/snpEff_summary.html \
        $infile \
        | bgzip -@ 10 -c - > \
        $outfile
    
    echo "Sim file annotated"

else
    echo "ERROR: Invalid species: $speies"
    exit 1
fi


echo "Finished annotating"
echo "Indexing..."

bcftools index $outfile

echo "Indexing complete"
