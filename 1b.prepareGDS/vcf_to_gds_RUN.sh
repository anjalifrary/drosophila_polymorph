#!/usr/bin/env bash
#SBATCH -J vcf2gds         # Job name
#SBATCH --ntasks=1         # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1               # Run on one node
#SBATCH -t 0-10:00         # 10 hours runtime
#SBATCH --mem=100G         # Memory per node
#SBATCH -o /scratch/ejy4bu/err_outs/GDS/vcf2gds.%A_%a.out  # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/GDS/vcf2gds.%A_%a.err  # Standard error
#SBATCH -p standard        # Partition
#SBATCH --account=berglandlab


module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1
export R_LIBS_USER=~/Rlibs

Rscript vcf_to_gds.R