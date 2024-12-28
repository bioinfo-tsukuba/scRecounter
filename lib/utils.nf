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