#!/usr/bin/env bash
#SBATCH -J gowinda         # Job name
#SBATCH --ntasks=1         # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1               # Run on one node
#SBATCH -t 0-10:00         # 10 hours runtime
#SBATCH --mem=70G         # Memory per node
#SBATCH -o /scratch/ejy4bu/err_outs/gowinda/gowinda.%A_%a.out  # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/gowinda/gowinda.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab
#SBATCH --array=0-10

mkdir -p /scratch/ejy4bu/err_outs/gowinda/

background=/scratch/ejy4bu/drosophila/gowinda/MAF5/new_6-29-26/background_all_snps.txt
# background=/scratch/ejy4bu/drosophila/gowinda/background_classed_snps.txt
# candidate_snp=/scratch/ejy4bu/drosophila/gowinda/candidate_snp_AB.txt
gtf_file=/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/dmel-all-r6.67.gtf
# go_file=/scratch/ejy4bu/drosophila/gowinda/flybase_go.txt
go_file=/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/flybase_gaf_go.txt
jar_file=/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/Gowinda-1.12.jar

SUFFICES=("A" "B" "AB" "FGOP" "XY" "FGOPXY" "ABFGOPXY")
SUFFIX=${SUFFICES[$SLURM_ARRAY_TASK_ID]}

echo "SUFFIX=$SUFFIX"
echo "--candidate-snp-file /scratch/ejy4bu/drosophila/gowinda/MAF5/new_6-29-26/candidate_snp_${SUFFIX}.txt"

# gene mode
# gene definition = gene (no exonic snps = not useful for updownstream)

java -Xmx8g -jar $jar_file \
  --snp-file $background \
  --candidate-snp-file /scratch/ejy4bu/drosophila/gowinda/MAF5/new_6-29-26/candidate_snp_${SUFFIX}.txt \
  --gene-set-file $go_file \
  --annotation-file $gtf_file \
  --simulations 1000000 \
  --gene-definition gene \
  --threads 10 \
  --mode snp \
  --output-file /scratch/ejy4bu/drosophila/gowinda/MAF5/new_6-29-26/results/gowinda_${SUFFIX}_snp_allBackground.txt



# once done exploring, set simulations to 1000000 (1M for final results)

# java -Xmx8g -jar /scratch/ejy4bu/drosophila/gowinda/Gowinda-1.12.jar \
#   --snp-file $background \
#   --candidate-snp-file $candidate_snp \
#   --gene-set-file $go_file \
#   --annotation-file $gtf_file \
#   --simulations 1000000 \
#   --gene-definition gene \
#   --threads 10 \
#   --mode snp \
#   --output-file /scratch/ejy4bu/drosophila/gowinda/results/final/gowinda_${suffix}_snp_allBackground.txt


  
# java -Xmx8g -jar /scratch/ejy4bu/drosophila/gowinda/Gowinda-1.12.jar \
#   --snp-file $background \
#   --candidate-snp-file $candidate_snp \
#   --gene-set-file $go_file \
#   --annotation-file $gtf_file \
#   --simulations 100000 \
#   --gene-definition updownstream2000 \
#   --threads 10 \
#   --mode snp \
#   --output-file /scratch/ejy4bu/drosophila/gowinda/results/gowinda_${suffix}_updown2k_classedBackground.txt

###### ARRAY run of all groups:
# #SBATCH --array=0-5


# suffix=$(basename "$candidate_snp" .txt)
# suffix=${suffix#candidate_snp_}

