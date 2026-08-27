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
#SBATCH --account berglandlab

set -euo pipefail

module load bcftools
# gff coord are 0-based, half-open = subtract 1 from the start pos, leave end pos as is 
# bed coord are 1-based, inclusive 

## sim
wm_dust="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/repeat/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna.wm.dust.bed"
rpt_mask="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/repeat/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna.out.gff"
filter_bed="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/repeat/dsim3.repeatMask_wmdust_combined.bed"
outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/"
vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.vcf.gz"
out_vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.vcf.gz"

## mel
rpt_mask="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/repeat/dmel-all-chromosome-r6.12.fasta.out.gff"
wm_dust="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/repeat/dmel-all-chromosome-r6.12.fasta.wm.dust.bed"
filter_bed="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/repeat/dm6.repeatMask_wmdust_combined.bed"
outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/"
vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.vcf.gz"
out_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.vcf.gz"


# r1_mel="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dmel-all-chromosome-r6.12.fasta.gz"
# r2_mel="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna"
### chr lengths are the same... so proceeding using these files 

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


# outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/"
# vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.vcf.gz"
# out_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.repeatmasked.wmdust.vcf.gz"
# bcftools view \
#     -T ^dmel.repeats_lowcomplexity.bed \
#     -Oz \
#     -o $out_vcf \
#     $vcf

# bcftools index -t $out_vcf


# cat /project/berglandlab/multispecies_endemism/data/repeat_masking_finalresult/dmel-all-chromosome-r6.12.fasta.wm.dust.bed | \
#     awk '{print($0"\twmDust")}' > /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.bed

#     cat /project/berglandlab/multispecies_endemism/data/repeat_masking_finalresult/dmel-all-chromosome-r6.12.fasta.out.gff | \
#     awk '{print($1"\t"$4"\t"$5"\trepasker") }' >> /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.bed
#     sed -i 's/-1/0/g' /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.bed

#     cat /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.bed | awk '{
#     if($2<$3) print $0
#     if($3<$2) print $1"\t"$3"\t"$2"\t"$4"flip"
#     }' > /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.flip.bed

