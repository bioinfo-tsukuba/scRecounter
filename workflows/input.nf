// Input workflow for processing paired-end reads
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

// Merge reads by sample; account for any differences in compression; check sequence formatting
process MERGE_READS {
    conda "envs/read_qc.yml"
    
    input:
    tuple val(sample), path("*_read1.fq"), path("*_read2.fq")

    output:
    tuple val(sample), path("${sample}_R1.fq"), path("${sample}_R2.fq")

    script:
    """
    seqkit seq *_read1.fq > ${sample}_R1.fq
    seqkit seq *_read2.fq > ${sample}_R2.fq
    """

    stub:
    """
    touch ${sample}_R1.fq ${sample}_R2.fq
    """
}
