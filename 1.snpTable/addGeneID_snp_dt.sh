#!/usr/bin/env bash
#SBATCH -J addGeneID       # Job name
#SBATCH --ntasks=1         # Single task per job
#SBATCH --cpus-per-task=1  # Number of CPU cores per task
#SBATCH -N 1               # Run on one node
#SBATCH -t 0-10:00         # 10 hours runtime
#SBATCH --mem=100G         # Memory per node
#SBATCH -o /scratch/ejy4bu/err_outs/geneID.%A_%a.out  # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/geneID.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab_standard


cd 1.snpTable/
export R_LIBS_USER="/sfs/gpfs/tardis/home/ejy4bu/R/goolf/4.5/"
module load gcc/11.4.0  openmpi/4.1.4 icu R/4.5.0
module load bcftools/1.23
Rscript addGeneID_snp_dt.R