#!/usr/bin/env bash
#SBATCH -J gowinda         # Job name
#SBATCH --ntasks=1         # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1               # Run on one node
#SBATCH -t 0-10:00         # 10 hours runtime
#SBATCH --mem=100G         # Memory per node
#SBATCH -o /scratch/ejy4bu/err_outs/gowinda.%A_%a.out  # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/gowinda.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab

total_snp=/scratch/ejy4bu/drosophila/gowinda/total_snp.txt
candidate_snp=/scratch/ejy4bu/drosophila/gowinda/candidate_snp_A.txt
gtf_file=/scratch/ejy4bu/drosophila/gowinda/dmel-all-r6.67.gtf
go_file=/scratch/ejy4bu/drosophila/gowinda/flybase_go.txt


java -Xmx8g -jar Gowinda.jar \
  --snp-file total_snp \
  --candidate-snp-file candidate_snp_A \
  --gene-set-file go_file \
  --annotation-file gtf_file \
  --simulations 1000000 \
  --gene-definition gene \
  --threads 16 \
  --mode gene \
  --output-file /scratch/ejy4bu/drosophila/gowinda/gowinda_A.txt