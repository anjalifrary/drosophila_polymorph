# ijob -A berglandlab -c24 -p standard --mem=64G

module load apptainer/1.3.4


cp /project/berglandlab/dylanhighland/ts_polymorphisms/drosophila_tsp/nflo_lastz/dsim_v3.1_to_dmel_v6/chainnet/liftover.chain \
/scratch/aob2x/dsim_v3.1_to_dmel_v6.chain
sed 's/NC_052520.2/sim_2L/g' -i /scratch/aob2x/dsim_v3.1_to_dmel_v6.chain
sed 's/NC_052521.2/sim_2R/g' -i /scratch/aob2x/dsim_v3.1_to_dmel_v6.chain
sed 's/NC_052522.2/sim_3L/g' -i /scratch/aob2x/dsim_v3.1_to_dmel_v6.chain
sed 's/NC_052523.2/sim_3R/g' -i /scratch/aob2x/dsim_v3.1_to_dmel_v6.chain
sed 's/NC_052524.2/sim_4/g' -i  /scratch/aob2x/dsim_v3.1_to_dmel_v6.chain
sed 's/NC_052525.2/sim_X/g' -i  /scratch/aob2x/dsim_v3.1_to_dmel_v6.chain


https://github.com/maurya-anand/liftover
apptainer build --docker-login /scratch/aob2x/bcftools_liftover.sif docker://ghanand/liftover:1.0.0
singularity shell /scratch/aob2x/bcftools_liftover.sif
bcftools +liftover --help


#head -n 1000 /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.vcf >  /scratch/aob2x/test.vcf
#nano /scratch/aob2x/test.vcf
#
#cp /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.vcf /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.orig.vcf
#sed -i 's/Number=\./Number=A/g' /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.vcf

bcftools +liftover \
-Ob -o /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf \
/scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.vcf -- \
-f /project/berglandlab/Dmel_genomic_resources/References/DESTv3_dmelholo/parts/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.cleanNames.fna \
-s /project/berglandlab/Dmel_genomic_resources/References/DESTv3_dmelholo/parts/GCF_016746395.2_Prin_Dsim_3.1_genomic.cleanNames.fna \
-c /scratch/aob2x/dsim_v3.1_to_dmel_v6.chain \
--drop-tags FORMAT/FREQ,FORMAT/AD \
--write-src



module load samtools bcftools gcc/11.4.0 openmpi/4.1.4 R/4.3.1

bcftools sort /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf |
sed 's/Number=A/Number=\./g' |
bgzip -c - > /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf.gz


Rscript --vanilla ~/DESTv3/snpCalling_dev/scatter_gather_annotate/vcf2gds.R \
/scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.vcf.gz

mv /scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds.gz \
/scratch/aob2x/20Nov2025_sim_dest3/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.gds