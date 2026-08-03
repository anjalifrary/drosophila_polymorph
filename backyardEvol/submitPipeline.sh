#!/usr/bin/env bash
#
#SBATCH -J make_bam_pipeline    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=1 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 0-18:00        # 10 hours runtime
#SBATCH --mem=5G        # Memory per node
#SBATCH -o /scratch/ejy4bu/err_outs/be/pipeline.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/be/pipeline.%A_%a.err # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab

set -euo pipefail


SAMPLE_LIST="/scratch/ejy4bu/backyardEvolution/metadata/testSamples.txt"
# SAMPLE_LIST="/scratch/ejy4bu/backyardEvolution/metadata/allSamples.txt"

trim_job=$(sbatch --parsable \
    --array=1-$(wc -l < ${SAMPLE_LIST})%30 \
    ~/drosophila_polymorph/backyardEvol/1.removeAdapters.sh)

map_job=$(sbatch --parsable \
    --dependency=afterok:${trim_job} \
    --array=1-$(wc -l < ${SAMPLE_LIST})%20 \
    ~/drosophila_polymorph/backyardEvol/2.map.sh)

dup_job=$(sbatch --parsable \
    --dependency=afterok:${map_job} \
    --array=1-$(wc -l < ${SAMPLE_LIST})%50 \
    ~/drosophila_polymorph/backyardEvol/3.markdup.sh)

clip_job=$(sbatch --parsable \
    --dependency=afterok:${dup_job} \
    --array=1-$(wc -l < ${SAMPLE_LIST})%50 \
    ~/drosophila_polymorph/backyardEvol/4.bamUtil.sh)

echo "Pipeline submitted:"
echo "trim: $trim_job"
echo "map:  $map_job"
echo "dup:  $dup_job"
echo "clip: $clip_job"