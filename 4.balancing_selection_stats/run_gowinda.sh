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
candidate_snp=/scratch/ejy4bu/drosophila/gowinda/candidate_snp_B.txt
gtf_file=/scratch/ejy4bu/drosophila/gowinda/dmel-all-r6.67.gtf
go_file=/scratch/ejy4bu/drosophila/gowinda/flybase_go.txt
# go_file=/scratch/ejy4bu/drosophila/gowinda/gene_association.fb

echo "total_snp=$total_snp"
echo "candidate_snp=$candidate_snp"
echo "gtf_file=$gtf_file"
echo "go_file=$go_file"

java -Xmx8g -jar /scratch/ejy4bu/drosophila/gowinda/Gowinda-1.12.jar \
  --snp-file $total_snp \
  --candidate-snp-file $candidate_snp \
  --gene-set-file $go_file \
  --annotation-file $gtf_file \
  --simulations 1000000 \
  --gene-definition gene \
  --threads 10 \
  --mode gene \
  --output-file /scratch/ejy4bu/drosophila/gowinda/results/gowinda_B_gene.txt


  
java -Xmx8g -jar /scratch/ejy4bu/drosophila/gowinda/Gowinda-1.12.jar \
  --snp-file $total_snp \
  --candidate-snp-file $candidate_snp \
  --gene-set-file $go_file \
  --annotation-file $gtf_file \
  --simulations 1000000 \
  --gene-definition updownstream2000 \
  --gene-definition gene \
  --threads 10 \
  --mode gene \
  --output-file /scratch/ejy4bu/drosophila/gowinda/results/gowinda_B_updown2k_gene.txt
