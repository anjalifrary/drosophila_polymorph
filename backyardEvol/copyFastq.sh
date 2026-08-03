#!/usr/bin/env bash
#
#SBATCH -J copy_fastq # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/ejy4bu/err_outs/be/copy.%A_%a.out  # Std out
#SBATCH -e /scratch/ejy4bu/err_outs/be/copy.%A_%a.out  # Std error
#SBATCH -p standard
#SBATCH --account berglandlab

### to Run
# sbatch --array=1-$( wc -l < /scratch/ejy4bu/backyardEvolution/metadata/sim_samps.csv )%10 ~/drosophila_polymorph/backyardEvol/copyFastq.sh
# sbatch --array=1-$( wc -l < /scratch/ejy4bu/backyardEvolution/metadata/mel_samps.csv )%10 ~/drosophila_polymorph/backyardEvol/copyFastq.sh

# cp /project/berglandlab/alan/be_flies/04.sampleMetadata/mel_samps.csv /scratch/ejy4bu/backyardEvolution/metadata/

module load gcc/11.4.0 sratoolkit/3.1.1 
# module load gcc/11.4.0 sratoolkit/3.1.1 aspera-connect/4.2.8

wd=/scratch/ejy4bu/backyardEvolution/
if [ ! -d $wd ]; then
  mkdir $wd
fi

# metadata="/scratch/ejy4bu/backyardEvolution/metadata/mel_samps.csv"
metadata="/scratch/ejy4bu/backyardEvolution/metadata/sim_samps.csv"


sampName=$(tail -n +2 "$metadata" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d",")
fastq1=$(tail -n +2 "$metadata" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d",")
fastq2=$(tail -n +2 "$metadata" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f3 -d",")
dir=$(basename "$fastq1" | sed 's/_WKDL.*//')

fastq1=$(dirname $fastq1)/${dir}/$(basename $fastq1)
fastq2=$(dirname $fastq2)/${dir}/$(basename $fastq2)

echo $sampName
echo $dir
echo $fastq1
echo $fastq2


if [ ! -d "${wd}/fastq/${sampName}" ]; then
  mkdir -p ${wd}/fastq/${sampName}
fi

cp $fastq1 "${wd}/fastq/${sampName}/${sampName}_1.fq.gz"
cp $fastq2 "${wd}/fastq/${sampName}/${sampName}_2.fq.gz"