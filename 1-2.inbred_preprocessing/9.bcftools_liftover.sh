# ijob -A berglandlab_standard -c2 -p standard --mem=64G

module load apptainer/1.5.0
module load samtools bcftools gcc/11.4.0 openmpi/4.1.4 R/4.5.0

singularity shell /scratch/ejy4bu/drosophila/liftover/bcftools_liftover.sif

outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/dsim3.signor/"
input_vcf_dsim3="${outdir}/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.vcf.gz"
vcf_dm6="${outdir}/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.merged.biallelic.dm6.vcf.gz"

ref_dsim3=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna
ref_dm6=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna
# chain_dsim3_dm6=/project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/dsim_v3.1_to_dmel_v6.chain
chain_dsim3_dm6=/project/berglandlab/anjali/drosophila_polymorphism/liftover/liftover.dsim3.dm6.chain

reject_vcf="${outdir}/dsim3.signor.dsim3TOdm6.liftover_rejected.vcf.gz"

bcftools +liftover \
  -Oz -o $vcf_dm6 \
  $input_vcf_dsim3 -- \
  -f $ref_dm6 \
  -s $ref_dsim3 \
  -c $chain_dsim3_dm6 \
  --reject $reject_vcf \
  --reject-type z \
  --drop-tags FORMAT/FREQ,FORMAT/AD \
  --write-src

bcftools sort "$vcf_dm6" -Oz -o "${vcf_dm6%.vcf.gz}.sorted.vcf.gz"

bcftools index -t "${vcf_dm6%.vcf.gz}.sorted.vcf.gz"

echo "Liftover, sorting, indexing are complete."


### Inspect rejected records:
tmp_in=$(mktemp)
tmp_rej=$(mktemp)

bcftools query -f '%CHROM\n' "$input_vcf_dsim3" |
    sort | uniq -c |
    awk '{print $2 "\t" $1}' |
    sort -k1,1 > "$tmp_in"

bcftools query -f '%CHROM\n' "$reject_vcf" |
    sort | uniq -c |
    awk '{print $2 "\t" $1}' |
    sort -k1,1 > "$tmp_rej"

{
    printf 'CHROM\tINPUT\tREJECTED\tREJECT_PCT\n'
    join -t $'\t' -a 1 -e 0 -o '0,1.2,2.2' "$tmp_in" "$tmp_rej" |
        awk -F '\t' 'BEGIN {OFS="\t"}
        $2 > 0 {
            printf "%s\t%d\t%d\t%.3f\n", $1, $2, $3, 100*$3/$2
        }'
} | {
    read -r header
    printf '%s\n' "$header"
    sort -t $'\t' -k4,4nr
}

rm "$tmp_in" "$tmp_rej"
