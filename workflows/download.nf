include { readAccessions; } from '../lib/download.groovy'
include { joinReads; addStats; } from '../lib/utils.nf'
include { SRA_STAT } from '../lib/utils.nf'

workflow DOWNLOAD_WF {
    take:
    ch_accessions

    main:
    // Run prefetch & fastq-dump
    ch_fqdump = FASTQ_DUMP(ch_accessions)

    /// Merge logs
    FQDUMP_LOG_MERGE(ch_fqdump.log.collect())
    
    // Join R1 and R2 channels, which will filter out empty R2 records
    ch_fastq = joinReads(ch_fqdump.R1, ch_fqdump.R2)

    emit:
    fastq = ch_fastq
}

process FQDUMP_LOG_MERGE {
    publishDir file(params.output_dir) / "logs", mode: "copy", overwrite: true
    label "download_env"

    input:
    path "*_log.csv"

    output:
    path "fq-dump_summary.csv"

    script:
    """
    csv-merge.py --outfile fq-dump_summary.csv *_log.csv
    """

    stub:
    """
    touch fq-dump_summary.csv 
    """
}

process FASTQ_DUMP {
    label "download_env"

    input:
    tuple val(sample), val(accession), val(metadata), val(sra_file_size_gb)

    output:
    tuple val(sample), val(accession), val(metadata), path("reads/read_1.fastq.gz"), emit: "R1"
    tuple val(sample), val(accession), val(metadata), path("reads/read_2.fastq.gz"), emit: "R2", optional: true
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
    touch reads/${accession}_1.fastq.gz reads/${accession}_2.fastq.gz
    """
}
