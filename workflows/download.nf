include { readAccessions; joinReads } from '../lib/download.groovy'

workflow DOWNLOAD_WF {
    main:
    // Load accessions from file
    ch_accessions = readAccessions(params.accessions)

    // Set up vdb-config with GCP credentials, if provided
    VDB_CONFIG()

    // Run sra-stat
    SRA_STAT(ch_accessions, VDB_CONFIG.out)

    // Run prefetch & fast(er)q-dump
    if ( params.max_spots > 0 ){
        // Subsample reads
        ch_fqdump = FASTQ_DUMP(ch_accessions)
    } else {
        // Download all reads
        PREFETCH(ch_accessions, VDB_CONFIG.out)
        PREFETCH_LOG_MERGE(PREFETCH.out.log.collect())
        ch_fqdump = FASTERQ_DUMP(PREFETCH.out.sra)
    }
    /// Merge logs
    FQDUMP_LOG_MERGE(ch_fqdump.log.collect())
    
    // Join R1 and R2 channels, which will filter out empty R2 records
    ch_fastq = joinReads(ch_fqdump.R1, ch_fqdump.R2)

    emit:
    fastq = ch_fastq
    sra_stat = SRA_STAT.out
}

process FQDUMP_LOG_MERGE {
    publishDir file(params.output_dir) / "logs", mode: "copy", overwrite: true
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-download/sc-recounter-download:0.1.0"
    conda "envs/download.yml"

    input:
    path "*_log.csv"

    output:
    path "fq-dump.csv"

    script:
    """
    csv-merge.py --outfile fq-dump.csv *_log.csv
    """

    stub:
    """
    touch fq-dump.csv 
    """
}

process FASTERQ_DUMP {
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-download/sc-recounter-download:0.1.0"
    conda "envs/download.yml"
    scratch { sra_file.size() < 200.GB ? "ram-disk" : false }
    memory { sra_file.size() < 200.GB ? (sra_file.size() / 1e9).GB * (task.attempt + 1) + 6.GB : 32.GB * task.attempt }
    time { sra_file.size() < 200.GB ? 24.h * task.attempt : 48.h + 12.h * task.attempt }
    cpus 8
    maxRetries 2

    input:
    tuple val(sample), val(accession), val(metadata), path(sra_file) 

    output:
    tuple val(sample), val(accession), val(metadata), path("reads/read_1.fastq"), emit: "R1"
    tuple val(sample), val(accession), val(metadata), path("reads/read_2.fastq"), emit: "R2", optional: true
    path "reads/fq-dump_log.csv", emit: "log"

    script:
    """
    fq-dump.py \\
      --sample ${sample} \\
      --threads ${task.cpus} \\
      --bufsize 10MB \\
      --curcache 50MB \\
      --mem 5GB \\
      --temp TMP_FILES \\
      --min-read-length ${params.min_read_len} \\
      --maxSpotId ${params.max_spots} \\
      --outdir reads \\
      ${sra_file}

    # remove the temporary files
    rm -rf TMP_FILES
    """

    stub:
    """
    mkdir -p reads
    touch reads/${accession}_1.fastq reads/${accession}_2.fastq
    """
}

process FASTQ_DUMP {
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-download/sc-recounter-download:0.1.0"
    conda "envs/download.yml"

    input:
    tuple val(sample), val(accession), val(metadata)

    output:
    tuple val(sample), val(accession), val(metadata), path("reads/read_1.fastq"), emit: "R1"
    tuple val(sample), val(accession), val(metadata), path("reads/read_2.fastq"), emit: "R2", optional: true
    path "reads/fq-dump_log.csv", emit: "log"

    script:
    """
    fq-dump.py \\
      --sample ${sample} \\
      --threads ${task.cpus} \\
      --bufsize 10MB \\
      --curcache 50MB \\
      --mem 5GB \\
      --temp TMP_FILES \\
      --min-read-length ${params.min_read_len} \\
      --maxSpotId ${params.max_spots} \\
      --outdir reads \\
      ${accession}

    # remove the temporary files
    rm -rf TMP_FILES
    """

    stub:
    """
    mkdir -p reads
    touch reads/${accession}_1.fastq reads/${accession}_2.fastq
    """
}

process PREFETCH_LOG_MERGE{
    publishDir file(params.output_dir) / "logs", mode: "copy", overwrite: true
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-download/sc-recounter-download:0.1.0"
    conda "envs/download.yml"

    input:
    path "*_log.csv"

    output:
    path "prefetch.csv"

    script:
    """
    csv-merge.py --outfile prefetch.csv *_log.csv
    """

    stub:
    """
    touch prefetch.csv 
    """
}

process PREFETCH {
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-download/sc-recounter-download:0.1.0"
    conda "envs/download.yml"
    label "process_long"

    input:
    tuple val(sample), val(accession), val(metadata)
    val vdb_config

    output:
    tuple val(sample), val(accession), val(metadata), path("prefetch_out/${accession}/${accession}.sra"), emit: sra, optional: true
    path "prefetch_out/prefetch_log.csv",  emit: log

    script:
    """
    export GCP_PROJECT_ID="${params.gcp_project_id}"
    export GCP_SQL_DB_TENANT="${params.db_tenant}"
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_user}"

    prefetch.py \\
      --sample ${sample} \\
      --max-size 5000 \\
      --outdir prefetch_out \\
      ${accession}
    """

    stub:
    """
    mkdir -p prefetch_out/${accession}
    touch prefetch_out/${accession}/${accession}.sra prefetch_out/prefetch-log.csv
    """
}

process SRA_STAT {
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-download/sc-recounter-download:0.1.0"
    conda "envs/download.yml"

    input:
    tuple val(sample), val(accession), val(metadata)
    val vdb_config

    output:
    tuple val(sample), val(accession), path("sra-stat.csv")

    script:
    """
    sra-stat.py ${accession}
    """

    stub:
    """
    touch sra-stat.csv
    """
}

process VDB_CONFIG {
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-download/sc-recounter-download:0.1.0"
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
