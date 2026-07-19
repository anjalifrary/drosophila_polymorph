#!/usr/bin/env bash
#
#SBATCH -J postprocess # A single job name for the array
#SBATCH --cpus-per-task=10
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/postpipeline.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/postpipeline.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### OLD VERSION _ UPDATE BEFORE RUN
### moved genomics dbi section to new script for optimization / memory constraints

set -euo pipefail

module load gcc/11.4.0
module load gatk
module load bcftools
module load htslib

ref="/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna"
outdir="/scratch/ejy4bu/drosophila/inbred/combined_vcf/"
mkdir -p ${outdir}
gvcf_dir="/scratch/ejy4bu/drosophila/inbred/fastq/PRJNA318623"

#  /project/berglandlab/alan/privatePolymorphisms/simulans/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.vcf.gz | \

chroms=("2L" "2R" "3L" "3R" "4" "X")
chr=${CHROMS[$SLURM_ARRAY_TASK_ID]}

JAVAMEM=80G
CPU=10

echo "Processing Chromosome: ${chr}"

# module load sratoolkit/3.1.1
# module load bwa
# module load samtools
# module load picard
# module load gatk
# module load fastqc

### 1. get gVCFs
# # find all vcfs
# find ${gvcf_dir} -name "*.g.vcf.gz" \
#     | sort > ${outdir}/gvcf.list

# echo "Found $(wc -l < ${outdir}/gvcf.list) gVCFs"
# if [ "$(wc -l < ${outdir}/gvcf.list)" -eq 0 ]; then
#     echo "ERROR: No gVCFs found in ${gvcf_dir}!" >&2
#     exit 1
# fi

###### sample map
sample_map="${outdir}/dsim_sample_map.txt"
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    echo "Finding gVCFs and creating sample map..."
    find ${gvcf_dir} -name "*.g.vcf.gz" | sort > ${outdir}/gvcf.list
    
    # Generate GATK-compliant sample map: [SampleName]\t[Path]
    awk '{n=$0; sub(".*/", "", n); sub(/\..*/, "", n); print n "\t" $0}' \
        ${outdir}/gvcf.list > ${sample_map}
        
    echo "Found $(wc -l < ${sample_map}) gVCFs."
fi

### from signor

# /home/cmb-07/sn1/ssignor/jre1.7.0/bin/java -Xmx20g -jar /home/cmb-07/sn1/ssignor/next-gen_example_sim/bin/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
# 	-T GenotypeGVCFs \
# 	-R $REF_PFX \
# 	-L X \
#   [the list of --variant ..]
#   -o .vcf

# gatk CombineGVCFs \
#     -R ${ref} \
#     $(while read gvcf; do echo "--variant $gvcf"; done < ${outdir}/gvcf.list) \
#     -O ${outdir}/dsim3.signor.combined.g.vcf.gz
# echo "Done combining gVCFs"

### Consolidate vcfs 

work="${outdir}/dsim_genomicsdb_${chr}"
tmp="/scratch/ejy4bu/tmp/temp_dsim_gdbi_${chr}"
# error handling:
if [ -d "$work" ]; then
    echo "Removing existing GenomicsDB workspace..."
    rm -rf $work
fi

mkdir -p $tmp

echo "Running GenomicsDBImport..."

gatk GenomicsDBImport \
    -R ${ref} \
    --genomicsdb-workspace-path $work \
    --tmp-dir $tmp \
    --batch-size 50 \
    --sample-name-map ${sample_map} \
    --reader-threads $CPU \
    -L ${chr}

rm -rf "$tmpdir"

### Joint Genotyping
gatk GenotypeGVCFs \
    -R ${ref} \
    -V gendb://${outdir}/dsim_genomicsdb \
    -O ${outdir}/dsim3.signor.combined.raw.vcf.gz

### Normalize vcfs
bcftools norm \
    -f ${ref} \
    -m -both \
    -Oz \
    -o ${outdir}/dsim3.signor.combined.norm.vcf.gz \
    ${outdir}/dsim3.signor.combined.raw.vcf.gz

bcftools index ${outdir}/dsim3.signor.combined.norm.vcf.gz

### Filter variants (used GATK best practices hard filters.. ?)
echo "Filtering..."
gatk VariantFiltration \
    -R ${ref} \
    -V ${outdir}/dsim3.signor.combined.norm.vcf.gz \
    -O ${outdir}/dsim3.signor.combined.filtered.vcf.gz \
    --filter-name "QD2" \
    --filter-expression "QD < 2.0" \
    --filter-name "FS60" \
    --filter-expression "FS > 60.0" \
    --filter-name "MQ40" \
    --filter-expression "MQ < 40.0" \
    --filter-name "MQRankSum125" \
    --filter-expression "MQRankSum < -12.5" \
    --filter-name "ReadPosRankSum8" \
    --filter-expression "ReadPosRankSum < -8.0"

# haplotype score?? deprecated in gatk 4+

### Signor's script
# /home/cmb-07/sn1/ssignor/jdk1.7.0_55/bin/java -Xmx2g -jar /home/cmb-07/sn1/ssignor/next-gen_example_sim/bin/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \
# 	-T VariantFiltration \
#     -R $REF_PFX \
#     --filterExpression "QD < 2.0" --filterName "QD2" \
#     --filterExpression "FS > 60.0" --filterName "FS60" \
#     --filterExpression "MQ < 40.0" --filterName "MQ40" \
#     --filterExpression "HaplotypeScore > 13.0" --filterName "HaplotypeScore13" \
#     --filterExpression "MQRankSum < -12.5" --filterName "MQRankSum125" \
#     --filterExpression "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum8" \
#     -o $OUT_PFX/$OUT_PFX.SNPs_filt.vcf \
#     -V $OUT_PFX/$OUT_PFX.SNPs.vcf;
# ckFileSz "$OUT_PFX/$OUT_PFX.SNPs_filt.vcf";

# ### subset for passed variants only -- SKIP (soft-filtering only)
# bcftools view \
#     -f PASS \
#     -Oz \
#     -o ${outdir}/dsim3.signor.combined.filtered.pass.vcf.gz \
#     ${outdir}/dsim3.signor.combined.filtered.vcf.gz

# bcftools index ${outdir}/dsim3.signor.combined.filtered.pass.vcf.gz


### SnpEff 
# note - if decide to run on passed samples only, modify in/out files
echo "Annotating with SnpEff..."
SNPEFF=/project/berglandlab/multispecies_endemism/snpEFF/v4.3t/snpEff

java -Xmx32G \
    -jar ${SNPEFF}/snpEff.jar ann -formatEff \
    -v Dsim_v3.1 \
    -stats ${outdir}/snpEff_summary.html \
    ${outdir}/dsim3.signor.combined.filtered.vcf.gz \
    | bgzip -@ 10 -c - > \
    ${outdir}/dsim3.signor.combined.filtered.ann.eff.vcf.gz

bcftools index ${outdir}/dsim3.signor.combined.filtered.ann.eff.vcf.gz

echo "Complete...!"