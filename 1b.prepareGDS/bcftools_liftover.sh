# ijob -A berglandlab -c2 -p standard --mem=64G

module load apptainer/1.5.0
module load samtools bcftools gcc/11.4.0 openmpi/4.1.4 R/4.5.0

# grep ">" $file_name | head -20

##########################################################################
### dsim2 -> dsim3.1 then dsim3.1 -> dm6

# chain file dsim2 -> dsim3.1
cp /project/berglandlab/dylanhighland/ts_polymorphisms/drosophila_tsp/nflo_lastz/dsim_v2_to_dsim_v3.1/chainnet/liftover.chain \
/project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/dsim_v2_to_dsim_v3.1.chain

chain_dsim2_3=/project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/dsim_v2_to_dsim_v3.1.chain

# chain file dsim3.1 -> dm6
cp /project/berglandlab/dylanhighland/ts_polymorphisms/drosophila_tsp/nflo_lastz/dsim_v3.1_to_dmel_v6/chainnet/liftover.chain \
/project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/dsim_v3.1_to_dmel_v6.chain

chain_dsim3_dm6=/project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/dsim_v3.1_to_dmel_v6.chain

# grep ">" $chain_dsim2_3 | head -20
# grep ">" $chain_dsim3_dm6 | head -20


# ### input vcf; from Signor, zenodo
# cd /project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/
# curl -L -o zenodo_sim.vcf.zip "https://zenodo.org/records/154261/files/simulans_multisamp_all_chr.vcf.zip?download=1"
# unzip zenodo_sim.vcf.zip

input_vcf_dsim2_orig=/project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/simulans_multisamp_all_chr.vcf
# cp /project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/simulans_multisamp_all_chr.vcf /project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/simulans_multisamp_all_chr_orig.vcf

# grep -v "^#" $input_vcf_dsim2 | cut -f1 | sort -u | head -20

### reference files 

# dsim2 reference 
cp /project/berglandlab/dylanhighland/ts_polymorphisms/drosophila_tsp/fastas/dsim-mod.fasta \
/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim-mod_v2.fasta

ref_dsim2=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim-mod_v2.fasta

# dsim 3.1 reference
cp /project/berglandlab/Dmel_genomic_resources/References/DESTv3_dmelholo/parts/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna \
/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna
# /project/berglandlab/dylanhighland/ts_polymorphisms/working_data/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.fasta

ref_dsim3=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna

# dm6 reference
cp /project/berglandlab/Dmel_genomic_resources/References/DESTv3_dmelholo/parts/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna \
/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna

ref_dm6=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna

# fix fasta headers
sed -i 's/^>Dsim_Scf_/>/' $ref_dsim2
samtools faidx $ref_dsim2 # make index .fai 

# mtDNA exist?
grep -i "mt" $ref_dsim2 | grep ">"

# fix chain naming: strip Dsim_Scf_ on chain1 source side, rename NC_* -> sim_* on chain1 target side and chain2 source side (both share the same six substitutions)
sed -i 's/Dsim_Scf_//g' $chain_dsim2_3
sed -i 's/NC_052520.2/sim_2L/g; s/NC_052521.2/sim_2R/g; s/NC_052522.2/sim_3L/g; s/NC_052523.2/sim_3R/g; s/NC_052524.2/sim_4/g; s/NC_052525.2/sim_X/g' $chain_dsim2_3
sed -i 's/NC_052520.2/sim_2L/g; s/NC_052521.2/sim_2R/g; s/NC_052522.2/sim_3L/g; s/NC_052523.2/sim_3R/g; s/NC_052524.2/sim_4/g; s/NC_052525.2/sim_X/g' $chain_dsim3_dm6


head -1 $chain_dsim2_3
head -1 $chain_dsim3_dm6
grep ">" $ref_dsim2 | head -6
grep ">" $ref_dsim3 | head -6

### reheader vcf... fix incorrect contig lengths (og vcf header lengths don't match dsim2 reference)
input_vcf_dsim2=/project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/zenodo_sim.reheadered.vcf

bcftools reheader -f ${ref_dsim2}.fai -o $input_vcf_dsim2 $input_vcf_dsim2_orig

# confirm header now matches reference lengths
bcftools view -h $input_vcf_dsim2 | grep "##contig" | head -10
cat ${ref_dsim2}.fai | head -10

# drop the single nt masked-position mismatch and any others
# bcftools norm --check-ref e -f $ref_dsim2 $input_vcf_dsim2 -o /dev/null
bcftools norm -f $ref_dsim2 -c s $input_vcf_dsim2 -Oz -o /project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/zenodo_sim.reheadered.fixed.vcf

# https://github.com/maurya-anand/liftover  --docker-login
# docker pull ghcr.io/maurya-anand/liftover

mkdir -p /scratch/ejy4bu/drosophila/liftover/

# don't rebuild apptainer
apptainer build /scratch/ejy4bu/drosophila/liftover/bcftools_liftover.sif docker://ghcr.io/maurya-anand/liftover
# apptainer build --docker-login /scratch/ejy4bu/drosophila/liftover/bcftools_liftover.sif docker://ghanand/liftover:1.0.0

singularity shell /scratch/ejy4bu/drosophila/liftover/bcftools_liftover.sif
# bcftools +liftover --help

#### dsim2 -> dsim3.1 
# have to reassign once in apptainer:
# input_vcf_dsim2=/project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/simulans_multisamp_all_chr.vcf
input_vcf_dsim2=/project/berglandlab/anjali/drosophila_polymorphism/data_files/vcfs/zenodo_sim.reheadered.vcf
ref_dsim3=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna
ref_dsim2=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/dsim-mod_v2.fasta
ref_dm6=/project/berglandlab/anjali/drosophila_polymorphism/data_files/fastas/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna
chain_dsim3_dm6=/project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/dsim_v3.1_to_dmel_v6.chain
chain_dsim2_3=/project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/dsim_v2_to_dsim_v3.1.chain

vcf_dsim3=/project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/simulans_multisamp_all_chr_dsim3.1.vcf

# singularity exec /scratch/ejy4bu/drosophila/liftover/bcftools_liftover.sif bash -c "
bcftools +liftover \
  -Ob -o $vcf_dsim3 \
  $input_vcf_dsim2 -- \
  -f $ref_dsim3 \
  -s $ref_dsim2 \
  -c $chain_dsim2_3 \
  --drop-tags FORMAT/FREQ,FORMAT/AD \
  --write-src
# "

# bcftools view $vcf_dsim3 -Ov -o /project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/simulans_multisamp_all_chr_dsim3.1.plain.vcf
### dsim3.1 -> dm6
vcf_dm6=/project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/simulans_multisamp_all_chr_dm6.vcf

singularity exec /scratch/ejy4bu/drosophila/liftover/bcftools_liftover.sif bash -c "
bcftools +liftover \
  -Ob -o $vcf_dm6 \
  $vcf_dsim3 -- \
  -f $ref_dm6 \
  -s $ref_dsim3 \
  -c $chain_dsim3_dm6 \
  --drop-tags FORMAT/FREQ,FORMAT/AD \
  --write-src
"

bcftools view -h $vcf_dm6 \
| sed 's/Number=A/Number=./g' > /project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/header.txt

bcftools view -H $vcf_dm6 \
> /project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/body.txt

cat /project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/header.txt /project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/body.txt \
| bcftools sort -Oz \
-o /project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/simulans_multisamp_all_chr_dm6.fixed.vcf.gz

bcftools index /project/berglandlab/anjali/drosophila_polymorphism/data_files/liftover/simulans_multisamp_all_chr_dm6.fixed.vcf.gz

##########################################################################################
# DEST liftover
##########################################################################################
# input file: /scratch/ejy4bu/drosophila/liftover/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.vcf
# dm6 reference: /project/berglandlab/Dmel_genomic_resources/References/DESTv3_dmelholo/parts/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna
# dsim3 reference: /project/berglandlab/Dmel_genomic_resources/References/DESTv3_dmelholo/parts/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna
# output: /scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf
# chain: /scratch/ejy4bu/drosophila/liftover/dsim_v3.1_to_dmel_v6.chain

# sed 's/NC_052520.2/sim_2L/g' -i /scratch/ejy4bu/drosophila/liftover/dsim_v3.1_to_dmel_v6.chain
# sed 's/NC_052521.2/sim_2R/g' -i /scratch/ejy4bu/drosophila/liftover/dsim_v3.1_to_dmel_v6.chain
# sed 's/NC_052522.2/sim_3L/g' -i /scratch/ejy4bu/drosophila/liftover/dsim_v3.1_to_dmel_v6.chain
# sed 's/NC_052523.2/sim_3R/g' -i /scratch/ejy4bu/drosophila/liftover/dsim_v3.1_to_dmel_v6.chain
# sed 's/NC_052524.2/sim_4/g' -i  /scratch/ejy4bu/drosophila/liftover/dsim_v3.1_to_dmel_v6.chain
# sed 's/NC_052525.2/sim_X/g' -i  /scratch/ejy4bu/drosophila/liftover/dsim_v3.1_to_dmel_v6.chain

#head -n 1000 /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.vcf >  /scratch/aob2x/test.vcf
#nano /scratch/aob2x/test.vcf
#
#cp /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.vcf /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.orig.vcf
#sed -i 's/Number=\./Number=A/g' /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.vcf

# cp /project/berglandlab/alan/privatePolymorphisms/simulans/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.vcf.gz \
# /scratch/ejy4bu/drosophila/liftover/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.vcf.gz

# gunzip -c /scratch/ejy4bu/drosophila/liftover/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.vcf.gz \
# > /scratch/ejy4bu/drosophila/liftover/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.vcf

sed -i 's/Number=\./Number=A/g' /scratch/ejy4bu/drosophila/liftover/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.vcf

# bcftools +liftover \
# -Ob -o /scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf \
# /scratch/ejy4bu/drosophila/liftover/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.NOREP.ann.eff.vcf -- \
# -f /project/berglandlab/Dmel_genomic_resources/References/DESTv3_dmelholo/parts/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna \
# -s /project/berglandlab/Dmel_genomic_resources/References/DESTv3_dmelholo/parts/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna \
# -c /scratch/ejy4bu/drosophila/liftover/dsim_v3.1_to_dmel_v6.chain \
# --drop-tags FORMAT/FREQ,FORMAT/AD \
# --write-src

# Output: Lines   total/swapped/reference added/rejected: 5209147/874371/238430/871051

cp /project/berglandlab/anjali/drosophila_polymorphism/dest/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf.gz \
/scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf.gz

gunzip -c /scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf.gz \
> /scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf

### create sorted and zipped files for conversion to gds 
module load samtools bcftools gcc/11.4.0 openmpi/4.1.4 R/4.3.1

bcftools view -h /scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf \
| sed 's/Number=A/Number=./g' > /scratch/ejy4bu/drosophila/liftover/header.txt

bcftools view -H /scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf \
> /scratch/ejy4bu/drosophila/liftover/body.txt

cat /scratch/ejy4bu/drosophila/liftover/header.txt /scratch/ejy4bu/drosophila/liftover/body.txt \
| bcftools sort -Oz \
-o /scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.fixed.vcf.gz

bcftools index /scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.fixed.vcf.gz


# ### ALan's version: create sorted and zipped files for conversion to gds 
# module load samtools bcftools gcc/11.4.0 openmpi/4.1.4 R/4.3.1
# bcftools sort /scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf |
# sed 's/Number=A/Number=\./g' |
# bgzip -c - > /scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf.gz


# gzip -c /scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf \
# > /scratch/ejy4bu/drosophila/liftover/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf.gz


# Rscript --vanilla ~/DESTv3/snpCalling_dev/scatter_gather_annotate/vcf2gds.R \
# /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf.gz

# mv /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds.gz \
# /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds

