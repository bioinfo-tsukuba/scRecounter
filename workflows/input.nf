workflow INPUT_WF{
    main:
    // load csv and extract accessions
    ch_samples = Channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true, sep: ",")
        .map { row ->
            def req_columns = ["sample", "fastq_1", "fastq_2"]
            def miss_columns = req_columns.findAll { !row.containsKey(it) }
            if (miss_columns) {
                error "Missing columns in the input CSV file: ${miss_columns}"
            }
            return [row.sample, file(row.fastq_1), file(row.fastq_2)]
        }.groupTuple()
        .map { sample, fastq_1, fastq_2 ->
            return [sample, fastq_1.flatten(), fastq_2.flatten()]
        }

    // merge reads by sample
    MERGE_READS(ch_samples)
    
    emit:
    fastq = MERGE_READS.out
}

process MERGE_READS {
    input:
    tuple val(sample), path(fastq_1), path(fastq_2)

    output:
    tuple val(sample), path("${sample}_R1.fq.gz"), path("${sample}_R2.fq.gz")

    script:
    """
    cat ${fastq_1} > ${sample}_R1.fq.gz
    cat ${fastq_2} > ${sample}_R2.fq.gz
    """
}