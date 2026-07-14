#!/usr/bin/env bash
#
#SBATCH -J download_SRA # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/ejy4bu/err_outs/prefetch.%A_%a.out  # Std out
#SBATCH -e /scratch/ejy4bu/err_outs/prefetch.%A_%a.out  # Std error
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc/11.4.0 sratoolkit/3.1.1 
# module load gcc/11.4.0 sratoolkit/3.1.1 aspera-connect/4.2.8

wd=/scratch/ejy4bu/drosophila/inbred
if [ ! -d $wd ]; then
  mkdir $wd
fi
### run as: sbatch --array=1-$( wc -l < ~/CompEvoBio_modules/utils/getSRA/sras_2025.csv )%20 ~/CompEvoBio_modules/utils/getSRA/downloadSRA.sh
### cat /scratch/aob2x/compBio/logs/prefetch.52222298_*.out | grep -B1 "do not"
### cat /scratch/aob2x/compBio/logs/prefetch.3259341_3.out


SLURM_ARRAY_TASK_ID=1

metadata="/project/berglandlab/anjali/drosophila_polymorphism/data_files/metadata/SraRunTable.csv"

sranum=$(tail -n +2 "$metadata" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d",")
proj=$(tail -n +2 "$metadata" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f6 -d",")
sampName=$(tail -n +2 "$metadata" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f40 -d",")

echo $sampName " / " $sranum " / " $proj

### Sz100  /  SRR3585384  /  PRJNA318623

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
  --outfile ${wd}/fastq/${proj}/${sranum} \
  -e 10 \
  -p \
  --temp /scratch/ejy4bu/tmp \
  ${wd}/sra/${sranum}.sra

  ls -lh ${wd}/fastq/${proj}/${sranum}*

fi

if [ -f "${wd}/fastq/${proj}/${sranum}_1.fastq" ]; then
  gzip ${wd}/fastq/${proj}/${sranum}_1.fastq
  gzip ${wd}/fastq/${proj}/${sranum}_2.fastq
fi

if [ -f "${wd}/fastq/${proj}/${sranum}" ]; then
  gzip -c ${wd}/fastq/${proj}/${sranum} > ${wd}/fastq/${proj}/${sranum}.fastq.gz
  rm ${wd}/fastq/${proj}/${sranum}
fi

#rm /scratch/aob2x/fastq/${sranum}.sra
# cat /home/aob2x/CompEvoBio_modules/data/runs.csv | nl | grep "SRR12463313"