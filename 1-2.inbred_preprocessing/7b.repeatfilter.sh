### from alan

cat /project/berglandlab/multispecies_endemism/data/repeat_masking_finalresult/dmel-all-chromosome-r6.12.fasta.wm.dust.bed | \
    awk '{print($0"\twmDust")}' > /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.bed

    cat /project/berglandlab/multispecies_endemism/data/repeat_masking_finalresult/dmel-all-chromosome-r6.12.fasta.out.gff | \
    awk '{print($1"\t"$4"\t"$5"\trepasker") }' >> /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.bed
    sed -i 's/-1/0/g' /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.bed

    cat /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.bed | awk '{
    if($2<$3) print $0
    if($3<$2) print $1"\t"$3"\t"$2"\t"$4"flip"
    }' > /project/berglandlab/alan/be_flies/00.refgenomes/dmel.reps.flip.bed