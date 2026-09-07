#!/usr/bin/env bash
#
#SBATCH -J primaryChr # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 50G
#SBATCH -o /scratch/ejy4bu/err_outs/dest/primarychr.%A.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/dest/primarychr.%A.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab



# ### DEST mel:
# outdir="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/filtering/"
in_vcf="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.eff.vcf.gz"
# gap_vcf="${outdir}/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.eff.snpgap10.vcf.gz"
# snp_vcf="${outdir}/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.eff.snpgap10.snpsOnly.vcf.gz"

# ### DEST sim:
# outdir="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/filtering/"
# in_vcf="/scratch/ejy4bu/drosophila/DEST_remake/vcfs/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.vcf.gz"
# gap_vcf="${outdir}/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.snpgap10.vcf.gz"
# snp_vcf="${outdir}/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.snpgap10.snpsOnly.vcf.gz"
# chr_names="${outdir}/renameChr.txt"
# out_vcf="${outdir}/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.primaryChrs.vcf.gz"


module purge
module load bcftools 

# primaryChr="2L,2R,3L,3R,4,X"
bcftools query -f '%CHROM\n' "$in_vcf" | sort -u

# bcftools annotate \
#     --rename-chrs $chr_names \
#     -Oz \
#     -o $out_vcf \
#     $in_vcf


# bcftools view -h "$in_vcf" | grep '^##contig'
# bcftools view \
#     -r "$primary" \
#     -Oz \
#     -o "$out_vcf" \
#     "$in_vcf"

# bcftools index -t "$out_vcf"