include { readAccessions; joinReads } from '../lib/download.groovy'

workflow DOWNLOAD_WF {
    main:
    // Load accessions from file
    ch_accessions = readAccessions(params.accessions)

    // Set up vdb-config with GCP credentials, if provided
    VDB_CONFIG()

    // Run sra-stat
    SRA_STAT(ch_accessions, VDB_CONFIG.out)
    /// Merge by sample
    //SRA_STAT_MERGE(SRA_STAT.out.groupTuple())

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
    //ch_fastq = joinReads(ch_fqdump.R1, ch_fqdump.R2) // merge by sample
    ch_fastq = ch_fqdump.R1.join(ch_fqdump.R2, by: [0,1])
 
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
    tuple val(sample), val(accession), path(sra_file) 

    output:
    tuple val(sample), val(accession), path("reads/${accession}_1.fastq"), emit: "R1"
    tuple val(sample), val(accession), path("reads/${accession}_2.fastq"), emit: "R2", optional: true
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

process FASTQ_DUMP {
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-download/sc-recounter-download:0.1.0"
    conda "envs/download.yml"

    input:
    tuple val(sample), val(accession)

    output:
    tuple val(sample), val(accession), path("reads/${accession}_1.fastq"), emit: "R1"
    tuple val(sample), val(accession), path("reads/${accession}_2.fastq"), emit: "R2", optional: true
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
    tuple val(sample), val(accession)
    val vdb_config

    output:
    tuple val(sample), val(accession), path("prefetch_out/${accession}/${accession}.sra"), emit: sra, optional: true
    path "prefetch_out/prefetch_log.csv",  emit: log

    script:
    """
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

process SRA_STAT_MERGE{
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-download/sc-recounter-download:0.1.0"
    conda "envs/download.yml"

    input:
    tuple val(sample), path("sra-stat_*.csv")

    output:
    tuple val(sample), path("sra-stat-merged.csv")

    script:
    """
    csv-merge.py --sample $sample --outfile sra-stat-merged.csv sra-stat_*.csv
    """

    stub:
    """
    touch sra-stat-merged.csv
    """
}

process SRA_STAT {
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-download/sc-recounter-download:0.1.0"
    conda "envs/download.yml"

    input:
    tuple val(sample), val(accession)
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
