#!/usr/bin/env bash
#
#SBATCH -J indel10bp # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/filt_indel.%A.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/filt_indel.%A.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

set -euo pipefail


module load bcftools
# module load bedtools

outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/"
in_vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.vcf.gz"
gap_vcf="${outdir}/dsim3.signor/combined.norm.gatkfilt.snpgap10.vcf.gz"
snp_vcf="${outdir}/dsim3.signor/combined.norm.gatkfilt.snpgap10.snpsOnly.vcf.gz"
# out_vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.indel10filt.snps.vcf.gz"
# ref="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"

bcftools filter \
    --SnpGap 10 \
    "$in_vcf" \
    -Oz \
    -o "$gap_vcf"

bcftools index -t "$gap_vcf"

bcftools view \
    -v snps \
    "$gap_vcf" \
    -Oz \
    -o "$snp_vcf"

bcftools index -t "$snp_vcf"


# genome="${outdir}/dsim3.chrLengths.txt"

# indel_file="${outdir}/dsim3.indels.bed"
# indel_10bp="${outdir}/dsim3.indelsplus10.bed"

# cut -f1,2 "${ref}.fai" > "$genome"

# # update later:
# # Get indels and convert to BED
# bcftools query -f '%CHROM\t%POS0\t%REF\t%ALT\n' "$in_vcf" |
# awk 'BEGIN{OFS="\t"} length($3) != length($4) {
#     print $1, $2, $2 + length($3)
# }' |
# bedtools sort -i - |
# bedtools slop -i - -g "$genome" -b 10 |
# bedtools merge -i - > "$indels10"

# # Remove variants within 10 bp of indels
# bcftools view \
#     -T ^"$indels10" \
#     -Oz \
#     -o "${outdir}/DGRP2.no_indel10.vcf.gz" \
#     "$in_vcf"

# bcftools index -t "${outdir}/DGRP2.no_indel10.vcf.gz"

# # Keep biallelic SNPs
# bcftools view \
#     -v snps \
#     -m2 -M2 \
#     -Oz \
#     -o "$out_vcf" \
#     "${outdir}/DGRP2.no_indel10.vcf.gz"

# bcftools index -t "$out_vcf"

# echo "Done:"
# echo "$out_vcf"





# bcftools view -v indels -H "$in_vcf" |
# awk 'BEGIN{OFS="\t"} {
#     print $1, $2-1, $2+length($4)-1, $4, $5
# }' \
# > $indel_file

# cut -f1,2 "${ref}.fai" > genome.txt

# bedtools slop \
#     -i $indel_file \
#     -g "${ref}.fai" \
#     -b 10 \
# > $indel_10bp

# bcftools query \
#     -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\n' \
#     "$in_vcf" |
# awk 'BEGIN{OFS="\t"} {
#     if (length($4) != length($5)) print $1,$2,$3,$4,$5
# }' \
# > $bed_file