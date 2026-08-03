




### 3. bwa mem single end mode -> to sorted bam
if [ ! -f "${SAMPLE_DIR}/bam/${samp_name}.sorted.bam" ]; then
    echo "mapping/sorting sample ${samp_name}. "
    bwa mem \
        -t 10 \
        -K 100000000 \
        -Y \
        -R "@RG\tID:${samp_name}\tSM:${samp_name}\tPL:ILLUMINA\tLB:${samp_name}" \
        ${ref} \
        ${SAMPLE_DIR}/fastp/${samp_name}.trimmed.fastq.gz \
        | samtools view -uh -q 20 \
        | samtools sort --threads 2 -o ${SAMPLE_DIR}/bam/${samp_name}.sorted.bam -

    samtools index ${SAMPLE_DIR}/bam/${samp_name}.sorted.bam
    echo "finished mapping $samp_name"
fi


bwa mem \
    -t 8 \
    reference.fa \
    sample_R1.trim.fastq.gz \
    sample_R2.trim.fastq.gz |
samtools view -b |
samtools sort -@8 -o sample.sorted.bam

samtools index sample.sorted.bam