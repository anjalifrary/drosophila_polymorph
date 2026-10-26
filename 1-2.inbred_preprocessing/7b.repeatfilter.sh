#!/usr/bin/env bash
#
#SBATCH -J repeatfilter # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/rpt_filt.%A.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/rpt_filt.%A.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

set -euo pipefail

module load bcftools
# bed coord are 0-based, half-open 
# gff coord are 1-based, inclusive = subtract 1 from the start pos, leave end pos as is 

## sim
outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/"

wm_dust="${outdir}/repeat/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna.wm.dust.bed"
rpt_mask="${outdir}/repeat/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna.out.gff"
filter_bed="${outdir}/repeat/dsim3.repeatMask_wmdust_combined.bed"

vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.vcf.gz"
out_vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.vcf.gz"

# ## mel
# outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/"

# rpt_mask="${outdir}/repeat/dmel-all-chromosome-r6.12.fasta.out.gff"
# wm_dust="${outdir}/repeat/dmel-all-chromosome-r6.12.fasta.wm.dust.bed"
# filter_bed="${outdir}/repeat/dm6.repeatMask_wmdust_combined.bed"

# vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.vcf.gz"
# out_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.vcf.gz"

{
    awk 'BEGIN{OFS="\t"} !/^#/ {
        print $1,$2,$3,"wmDust"
    }' "$wm_dust"

    awk 'BEGIN{OFS="\t"} !/^#/ {
        print $1,$4-1,$5,"RepeatMasker"
    }' "$rpt_mask"
} |
sort -k1,1 -k2,2n > $filter_bed
# sort by chr then by pos (n=numeric)


echo "created bed for filtering"



echo "filtering vcf"
bcftools view \
    -T ^$filter_bed \
    -Oz \
    -o $out_vcf \
    $vcf

bcftools index -t $out_vcf
echo "completed filtering & indexing"
