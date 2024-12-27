include { readAccessions; joinReads } from '../lib/download.groovy'

workflow DOWNLOAD_WF {
    take:
    ch_accessions

    main:
    // read accessions
    ch_accessions = readAccessions(ch_accessions)

    // Run sra-stat
    SRA_STAT(ch_accessions)

    // Run prefetch & fast(er)q-dump
    if ( params.max_spots > 0 ){
        // Subsample reads
        ch_fqdump = FASTQ_DUMP(ch_accessions)
    } else {
        // Download all reads
        PREFETCH(ch_accessions)
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
    label "download_env"

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
    label "download_env"
    memory { 
        def baseMem = 16.GB
        def additionalMem = sra_file.size() / 10e9
        baseMem + (additionalMem * task.attempt).GB
    }
    time { 
        def baseTime = 4.h
        def additionalHours = sra_file.size() / 20e9
        baseTime + (additionalHours * task.attempt).h
    }
    disk { 
        def baseSize = 100.GB
        def additionalSize = sra_file.size() * (5 + 5 * task.attempt)
        additionalSize < baseSize ? baseSize : additionalSize
    }
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
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    fq-dump.py \\
      --sample ${sample} \\
      --accession ${accession} \\
      --threads ${task.cpus} \\
      --bufsize 10MB \\
      --curcache 50MB \\
      --mem 5GB \\
      --temp TMP_FILES \\
      --min-read-length ${params.min_read_len} \\
      --maxSpotId ${params.max_spots} \\
      --outdir reads \\
      ${sra_file}
    """

    stub:
    """
    mkdir -p reads
    touch reads/${accession}_1.fastq reads/${accession}_2.fastq
    """
}

process FASTQ_DUMP {
    label "download_env"

    input:
    tuple val(sample), val(accession), val(metadata)

    output:
    tuple val(sample), val(accession), val(metadata), path("reads/read_1.fastq"), emit: "R1"
    tuple val(sample), val(accession), val(metadata), path("reads/read_2.fastq"), emit: "R2", optional: true
    path "reads/fq-dump_log.csv", emit: "log"

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    fq-dump.py \\
      --sample ${sample} \\
      --accession ${accession} \\
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
    label "download_env"

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
    label "download_env"
    label "process_long"

    input:
    tuple val(sample), val(accession), val(metadata)

    output:
    tuple val(sample), val(accession), val(metadata), path("prefetch_out/${accession}/${accession}.{sra,sralite,sharq}"), emit: sra, optional: true
    path "prefetch_out/prefetch_log.csv",  emit: log

    script:
    def gcp_download = params.executor == "google-batch" ? "--gcp-download" : ""
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    prefetch.py ${gcp_download} \\
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
    label "download_env"

    input:
    tuple val(sample), val(accession), val(metadata)

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

/*
process VDB_CONFIG {
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-download/sc-recounter-download:0.1.0"
    conda "envs/download.yml"

    output:
    path "config.mkfg"

    script:
    """
    if [[ "${params.executor}" == "google-batch" ]]; then
        vdb-config --root \\
          --accept-gcp-charges \\
          > config.mkfg
    else
        vdb-config > config.mkfg
    fi
    """
}
*/

/*
vdb-config --report-cloud-identity yes


        vdb-config --root \
          -s /repository/remote/main/GCP/public/root="gs://sra-pub-caching" \
          -s /repository/remote/protected/GCP/root="gs://sra-pub-caching" 


 vdb-config --root --accept-gcp-charges yes > config.mkfg

 */