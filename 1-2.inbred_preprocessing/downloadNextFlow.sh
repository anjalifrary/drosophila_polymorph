# https://nf-co.re/sarek/3.9.0/

java -version
# if older than java 17 do: 
## module load java/17 

cd $HOME/.local/bin

# download and install
curl -s https://get.nextflow.io | bash

# make binary executable
chmod +x nextflow

./nextflow info

echo 'export PATH=$HOME/.local/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

which nextflow
nextflow info

### to update
# nextflow self-update


### to run
module load java/17
module load apptainer
# export NXF_APPTAINER_CACHEDIR=/scratch/ejy4bu/tmp/apptainer_cache
# export APPTAINER_TMPDIR=/scratch/ejy4bu/tmp/apptainer_tmp
# export APPTAINER_CACHEDIR=/scratch/ejy4bu/tmp/apptainer_cache

### demo
nextflow run nf-core/demo -profile test,apptainer --outdir results
ls results/

# input file (CSV)
    # Column	Description
    # sample	Unique sample identifier
    # fastq_1	Path to gzipped first-read FastQ file
    # fastq_2	Path to gzipped second-read FastQ file (leave empty for single-end)

nextflow run nf-core/demo \
  -profile singularity \
  --input samplesheet.csv \
  --outdir results