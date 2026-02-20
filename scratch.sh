ijob -A berglandlab -c2 -p standard --mem=40G
export R_LIBS_USER="/sfs/gpfs/tardis/home/ejy4bu/R/goolf/4.5/"
module load gcc/11.4.0  openmpi/4.1.4 icu R/4.5.0
module load R/4.5.0
R


/// 
# !/usr/bin/env bash
# SBATCH -J jobname         # Job name
# SBATCH --ntasks=1         # Single task per job
# SBATCH --cpus-per-task=10 # Number of CPU cores per task
# SBATCH -N 1               # Run on one node
# SBATCH -t 0-10:00         # 10 hours runtime
# SBATCH --mem=100G         # Memory per node
# SBATCH -o /scratch/ejy4bu/err_outs/name.%A_%a.out  # Standard output
# SBATCH -e /scratch/ejy4bu/err_outs/name.%A_%a.err  # Standard error
# SBATCH -p standard       # Partition
# SBATCH --account=berglandlab