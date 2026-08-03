



### 4. dedup mark duplicates (GATK)

if [ ! -f "${SAMPLE_DIR}/bam/${samp_name}.markdup.bam" ]; then
    echo "dedup sample ${samp_name}. "
    java -Xmx45G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        I=${SAMPLE_DIR}/bam/${samp_name}.sorted.bam  \
        O=${SAMPLE_DIR}/bam/${samp_name}.markdup.bam \
        M=${SAMPLE_DIR}/bam/${samp_name}.markdup.metrics.txt \
        CREATE_INDEX=true
fi
# should i remove duplicates or just flag?


picard MarkDuplicates \
    INPUT=sample.sorted.bam \
    OUTPUT=sample.markdup.bam \
    METRICS_FILE=sample.metrics.txt \
    CREATE_INDEX=true