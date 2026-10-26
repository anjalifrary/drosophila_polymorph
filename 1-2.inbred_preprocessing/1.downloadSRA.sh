#!/usr/bin/env bash
#
#SBATCH -J download_SRA # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/ejy4bu/err_outs/SRA/prefetch.%A_%a.out  # Std out
#SBATCH -e /scratch/ejy4bu/err_outs/SRA/prefetch.%A_%a.out  # Std error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

### to Run
# sbatch --array=1-$( wc -l < /project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/SraRunTable.csv )%10 ~/drosophila_polymorph/1-2.inbred_dataProcessing/1.downloadSRA.sh

module load gcc/11.4.0 sratoolkit/3.1.1 

wd=/scratch/ejy4bu/drosophila/inbred
if [ ! -d $wd ]; then
  mkdir $wd
fi

metadata="/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/SraRunTable.csv"

sranum=$(tail -n +2 "$metadata" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d",")
proj=$(tail -n +2 "$metadata" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f6 -d",")
sampName=$(tail -n +2 "$metadata" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f40 -d",")

echo $sampName " / " $sranum " / " $proj

if [ ! -d "${wd}/fastq/${proj}" ]; then
  mkdir -p ${wd}/fastq/${proj}
fi


if ls ${wd}/fastq/${proj}/${sranum}*fastq.gz 1> /dev/null 2>&1; then
    echo "files do exist"
else
  echo "files do not exist"

  echo "force re-download"
  prefetch \
  -o ${wd}/sra/${sranum}.sra \
  -p \
  ${sranum}

  fasterq-dump \
  --split-files \
  --split-3 \
  -O ${wd}/fastq/${proj}/${sranum} \
  -e 10 \
  -p \
  --temp /scratch/ejy4bu/tmp \
  ${wd}/sra/${sranum}.sra

  ls -lh ${wd}/fastq/${proj}/${sranum}*

fi
