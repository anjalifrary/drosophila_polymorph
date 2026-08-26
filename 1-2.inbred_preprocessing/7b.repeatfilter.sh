### from alan

wm_dust="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/repeat/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna.wm.dust.bed"
rpt_mask="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/repeat/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna.out.gff"

{
    awk 'BEGIN{OFS="\t"} !/^#/ {
        print $1,$2,$3,"wmDust"
    }' "$wm_dust"

    awk 'BEGIN{OFS="\t"} !/^#/ {
        print $1,$4-1,$5,"RepeatMasker"
    }' "$rpt_mask"
} |
sort -k1,1 -k2,2n |
# bedtools merge -i - \
#     > dsim.repeats_lowcomplexity.bed



outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/"
vcf="${outdir}/dsim3.signor.combined.norm.gatkfilt.vcf.gz"
out_vcf="${outdir}//dsim3.signor.combined.norm.gatkfilt.repeatmasked.wmdust.vcf.gz"

bcftools view \
    -T ^dsim.repeats_lowcomplexity.bed \
    -Oz \
    -o $out_vcf \
    $vcf

bcftools index -t $out_vcf


# outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/DGRP2/"
# vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.vcf.gz"
# out_vcf="${outdir}/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.repeatmasked.wmdust.vcf.gz"
# bcftools view \
#     -T ^dmel.repeats_lowcomplexity.bed \
#     -Oz \
#     -o $out_vcf \
#     $vcf

# bcftools index -t $out_vcf


cat /project/berglandlab/multispecies_endemism/data/repeat_masking_finalresult/dmel-all-chromosome-r6.12.fasta.wm.dust.bed | \
    awk '{print($0"\twmDust")}' > /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.bed

    cat /project/berglandlab/multispecies_endemism/data/repeat_masking_finalresult/dmel-all-chromosome-r6.12.fasta.out.gff | \
    awk '{print($1"\t"$4"\t"$5"\trepasker") }' >> /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.bed
    sed -i 's/-1/0/g' /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.bed

    cat /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.bed | awk '{
    if($2<$3) print $0
    if($3<$2) print $1"\t"$3"\t"$2"\t"$4"flip"
    }' > /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.flip.bed

