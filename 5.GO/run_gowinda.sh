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
#SBATCH --array=0-7

mkdir -p /scratch/ejy4bu/err_outs/gowinda/

BG_dir=/scratch/ejy4bu/drosophila/gowinda/backgroundFiles/
CANDIDATE_root=/scratch/ejy4bu/drosophila/gowinda/candidateFiles/
OUTPUT_root=/scratch/ejy4bu/drosophila/gowinda/results

CANDIDATE_dir=MAF25filter_polyAF
CANDIDATE_id=25_polyAF

#maf_inputs <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.49)
# running 0.02, 0.01, 0.005, 5, 10, 15, 

gtf_file=/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/dmel-all-r6.67.gtf
go_file=/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/flybase_gaf_go.txt
jar_file=/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/Gowinda-1.12.jar

ls /scratch/ejy4bu/drosophila/gowinda/backgroundFiles/bg_sharedOnly_noMAF.txt


SUFFICES=("AB" "XY" "FGOPXY" "ABFGOPXY")
BACKGROUNDS=("bg_speciesSpecific_noMAF" "bg_sharedOnly_noMAF")

BG_INDEX=$((SLURM_ARRAY_TASK_ID / 4))
SUFFIX_INDEX=$((SLURM_ARRAY_TASK_ID % 4))

BACKGROUND=${BACKGROUNDS[$BG_INDEX]}
SUFFIX=${SUFFICES[$SUFFIX_INDEX]}

# SUFFIX=${SUFFICES[$SLURM_ARRAY_TASK_ID]}

echo "SUFFIX=$SUFFIX"
echo "BACKGROUND=$BACKGROUND"

background="${BG_dir}/${BACKGROUND}.txt"
candidate="${CANDIDATE_root}/${CANDIDATE_dir}/candidate_chrpos_${SUFFIX}_${CANDIDATE_id}.txt"
output="${OUTPUT_root}/${CANDIDATE_dir}/${BACKGROUND}/gowinda_${SUFFIX}_${CANDIDATE_id}.txt"

mkdir -p "${OUTPUT_root}/${CANDIDATE_dir}/${BACKGROUND}"

echo "bg=$background"
echo "candidate file=$candidate"
echo "output file=$output"

# gene mode (more conservative)
# gene definition = gene (no exonic snps = not useful for updownstream)

java -Xmx8g -jar $jar_file \
  --snp-file $background \
  --candidate-snp-file $candidate \
  --gene-set-file $go_file \
  --annotation-file $gtf_file \
  --simulations 1000000 \
  --gene-definition gene \
  --threads 10 \
  --mode gene \
  --output-file $output



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


# # background=/scratch/ejy4bu/drosophila/gowinda/MAF5/new_6-29-26/background_all_snps.txt
# background=/scratch/ejy4bu/drosophila/gowinda/background_all_snps.txt
# # background=/scratch/ejy4bu/drosophila/gowinda/background_classed_snps.txt
# # candidate_snp=/scratch/ejy4bu/drosophila/gowinda/candidate_snp_AB.txt
# gtf_file=/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/dmel-all-r6.67.gtf
# # go_file=/scratch/ejy4bu/drosophila/gowinda/flybase_go.txt
# go_file=/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/flybase_gaf_go.txt
# jar_file=/project/berglandlab/anjali/drosophila_polymorphism/gene_ontology/gowinda/Gowinda-1.12.jar
