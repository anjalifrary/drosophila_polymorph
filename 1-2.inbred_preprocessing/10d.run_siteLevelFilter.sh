#!/usr/bin/env bash
#SBATCH -J siteFilt_RD     # Job name
#SBATCH --ntasks=1         # Single task per job
#SBATCH --cpus-per-task=4 # Number of CPU cores per task
#SBATCH -N 1               # Run on one node
#SBATCH -t 0-10:00         # 10 hours runtime
#SBATCH --mem=50G         # Memory per node
#SBATCH -o /scratch/ejy4bu/err_outs/GDS/sampleFilt.%A_%a.out  # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/GDS/sampleFil.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab
#SBATCH --array=1-4

# NVAR=3949956 # mel
NVAR=2948855 # sim
# NVAR=1000 # test

CHUNK=$(( (NVAR + 3) / 4 ))

START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * CHUNK ))

# Don't go past the final variant
if [ $END -gt $NVAR ]; then
    END=$NVAR
fi

echo "Job ${SLURM_ARRAY_TASK_ID}"
echo "START = ${START}"
echo "END   = ${END}"


cd 1-2.inbred_preprocessing/
export R_LIBS_USER="/sfs/gpfs/tardis/home/ejy4bu/R/goolf/4.5/"
module load gcc/11.4.0  openmpi/4.1.4 icu R/4.5.0
Rscript 10d.siteLevelFilter.R $START $END



# junk:
# OUTDIR="/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter"

# NVAR=3949956
# CHUNK=$(( (NVAR + 3) / 4 ))

# head -n 1 "${OUTDIR}/site_RD_1_${CHUNK}.csv" \
#     > "${OUTDIR}/sim_site_RD.csv"

# for f in \
#     "${OUTDIR}/sim_site_RD_1_${CHUNK}.csv" \
#     "${OUTDIR}/sim_site_RD_$((CHUNK+1))_$((2*CHUNK)).csv" \
#     "${OUTDIR}/sim_site_RD_$((2*CHUNK+1))_$((3*CHUNK)).csv" \
#     "${OUTDIR}/sim_site_RD_$((3*CHUNK+1))_${NVAR}.csv"
# do
#     tail -n +2 "$f" >> "${OUTDIR}/sim_site_RD.csv"
# done