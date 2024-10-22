workflow DOWNLOAD_WF {
    main:
    // Load accessions from file
    ch_accessions = Channel
        .fromPath(params.accessions, checkIfExists: true)
        .splitCsv(header: true, sep: ",")
        .map { row ->
            def req_columns = ["sample", "accession"]
            def miss_columns = req_columns.findAll { !row.containsKey(it) }
            if (miss_columns) {
                error "Missing columns in the input CSV file: ${miss_columns}"
            }
            return [row.sample, row.accession]
        }

    // prefetch
    PREFETCH(ch_accessions)
    
    // fasterq-dump
    FASTERQ_DUMP(PREFETCH.out)

    // join R1 and R2 channels, which will filter out empty R2 records
    ch_fastq = FASTERQ_DUMP.out.R1.join(FASTERQ_DUMP.out.R2, by: [0, 1], remainder: true)
        .filter { sample, accession, r1, r2 -> 
            if(r2 == null) {
                println "Warning: Read 2 is empty for ${sample}-${accession}; skipping"
            }
            return r2 != null
        }
        .groupTuple()
        .map { sample, accession, fastq_1, fastq_2 ->
            return [sample, fastq_1.flatten(), fastq_2.flatten()]
        }

    // Merge by sample
    MERGE_READS(ch_fastq)

    emit:
    fastq = MERGE_READS.out
}

process MERGE_READS {
    conda "envs/read_qc.yml"
    scratch true
    
    input:
    tuple val(sample), path("*_read1.fq"), path("*_read2.fq")

    output:
    tuple val(sample), path("${sample}_R1.fq"), path("${sample}_R2.fq")

    script:
    """
    seqkit seq *_read1.fq > ${sample}_R1.fq
    seqkit seq *_read2.fq > ${sample}_R2.fq

    # remove the input files to save space
    if [[ "${params.keep_temp}"  != "true" ]]; then
        rm -f \$(readlink *_read1.fq) \$(readlink *_read2.fq)
    fi  
    """

    stub:
    """
    touch ${sample}_R1.fq ${sample}_R2.fq
    """
}

process FASTERQ_DUMP {
    conda "envs/download.yml"
    label "process_low"
    scratch true

    input:
    tuple val(sample), val(accession), path(sra_file)

    output:
    tuple val(sample), val(accession), path("reads/${accession}_1.fastq"), emit: "R1"
    tuple val(sample), val(accession), path("reads/${accession}_2.fastq"), emit: "R2", optional: true

    script:
    """
    f-dump.py \\
      --threads ${task.cpus} \\
      --bufsize 10MB \\
      --curcache 50MB \\
      --mem 5GB \\
      --temp TMP_FILES \\
      --maxSpotId ${params.max_spots} \\
      --outdir reads \\
      ${sra_file}

    # remove the input sra file
    if [[ "${params.keep_temp}"  != "true" ]]; then
        rm -f \$(readlink ${sra_file})
    fi
    """

    stub:
    """
    mkdir -p reads
    touch reads/${accession}_1.fastq reads/${accession}_2.fastq
    """
}

process PREFETCH {
    conda "envs/download.yml"
    label "process_low"
    scratch true

    input:
    tuple val(sample), val(accession)

    output:
    tuple val(sample), val(accession), path("prefetch_out/${accession}/${accession}.sra")

    script:
    """
    prefetch.py \\
      --max-size 5000 \\
      --outdir prefetch_out \\
      ${accession}
    """

    stub:
    """
    mkdir -p prefetch_out/${accession}
    touch prefetch_out/${accession}/${accession}.sra
    """
}

process VDB_DUMP_MERGE {
    publishDir file(params.outdir) / "download", mode: "copy", overwrite: true
    conda "envs/download.yml"
    label "process_low"

    input:
    path info

    output:
    path "info.csv"

    script:
    """
    vdb-dump_merge.py $info > info.csv
    """
}

process VDB_DUMP_INFO {
    conda "envs/download.yml"
    label "process_low"

    input:
    tuple val(sample), val(accession)

    output:
    tuple val(sample), val(accession), path("${sample}_${accession}_info.txt")

    script:
    """
    exit 1
    vdb-dump --info ${accession} > ${accession}_info.txt
    """
}

