process SRA_STAT {
    label "download_env"
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    disk 10.GB

    input:
    tuple val(sample), val(accession), val(download_url), val(metadata)

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