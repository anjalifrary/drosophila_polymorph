
module load apptainer/1.5.0
module load samtools bcftools gcc/11.4.0 openmpi/4.1.4 R/4.5.0

singularity shell /scratch/ejy4bu/drosophila/liftover/bcftools_liftover.sif


outdir="/project/berglandlab/alan/be_flies/07.vcf/"
input_vcf_dsim3="${outdir}/sim_be.r6.12.norm.filtered.repeatmasked.eff.biallelicSNP.fixedPloidy.vcf.gz"

ref_dsim3=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna
ref_dm6=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna
chain_dsim3_dm6=/project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/dsim_v3.1_to_dmel_v6.chain

# ouput vcf:
vcf_dm6="/scratch/ejy4bu/backyardEvolution/liftedOver/sim_be.r6.12.norm.filtered.repeatmasked.eff.biallelicSNP.fixedPloidy.dm6.vcf.gz"

bcftools +liftover \
  -Oz -o $vcf_dm6 \
  $input_vcf_dsim3 -- \
  -f $ref_dm6 \
  -s $ref_dsim3 \
  -c $chain_dsim3_dm6 \
  --drop-tags FORMAT/FREQ,FORMAT/AD \
  --write-src

bcftools sort "$vcf_dm6" -Oz -o "${vcf_dm6%.vcf.gz}.sorted.vcf.gz"

bcftools index -t "${vcf_dm6%.vcf.gz}.sorted.vcf.gz"
