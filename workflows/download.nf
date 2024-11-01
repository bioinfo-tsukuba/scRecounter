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
            // remove special characters from the sample name
            row.sample = row.sample.replaceAll("\\s", "_")
            return [row.sample, row.accession]
        }

    // Set up vdb-config with GCP credentials, if provided
    VDB_CONFIG()

    // Run prefetch
    PREFETCH(ch_accessions, VDB_CONFIG.out)
    
    // Run fast(er)q-dump
    FQ_DUMP(PREFETCH.out)

    // Join R1 and R2 channels, which will filter out empty R2 records
    ch_fastq = FQ_DUMP.out.R1.join(FQ_DUMP.out.R2, by: [0, 1], remainder: true)
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

    emit:
    fastq = ch_fastq
}

process FQ_DUMP {
    conda "envs/download.yml"
    scratch { sra_file.size() < 200.GB ? "ram-disk" : false }
    memory { sra_file.size() < 200.GB ? (sra_file.size() / 1e9).GB * (task.attempt + 1) + 6.GB : 32.GB * task.attempt }
    time { sra_file.size() < 200.GB ? 24.h * task.attempt : 48.h + 12.h * task.attempt }
    cpus 8
    maxRetries 2

    input:
    tuple val(sample), val(accession), path(sra_file) 

    output:
    tuple val(sample), val(accession), path("reads/${accession}_1.fastq"), emit: "R1"
    tuple val(sample), val(accession), path("reads/${accession}_2.fastq"), emit: "R2", optional: true

    script:
    """
    fq-dump.py \\
      --threads ${task.cpus} \\
      --bufsize 10MB \\
      --curcache 50MB \\
      --mem 5GB \\
      --temp TMP_FILES \\
      --maxSpotId ${params.max_spots} \\
      --outdir reads \\
      ${sra_file}

    # remove the temporary files
    rm -rf TMP_FILES

    # remove the input sra file
    if [[ "${params.keep_temp}" != "true" ]]; then
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
    label "process_long"

    input:
    tuple val(sample), val(accession)
    val vdb_config

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

process VDB_CONFIG {
    conda "envs/download.yml"

    output:
    val true

    script:
    """
    if [[ -f "${params.gcp_json}" ]] && [[ ! -z "${params.gcp_json}" ]]; then
        vdb-config set storage.gcs.service-account-file ${params.gcp_json}
    fi
    """
}
