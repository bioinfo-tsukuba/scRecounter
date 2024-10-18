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
            return [row.sample, row.fastq_1, row.fastq_2]
        }.groupTuple()
        .map { sample, fastq_1, fastq_2 ->
            return [sample, fastq_1.collect(), fastq_2.collect()]
        }

    emit:
    samples = ch_samples
}