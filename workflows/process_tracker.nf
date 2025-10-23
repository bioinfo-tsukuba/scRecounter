// Process tracker processes
process PROCESS_TRACKER_START {
    label "tracker_env"
    
    input:
    tuple val(sample), val(accession), val(download_url), val(metadata)
    val process_type
    val process_id

    script:
    def experiment_id = "${accession}"
    def organism = metadata.organism ?: ""
    """
    python3 ${projectDir}/bin/process_tracker_start.py \\
        --experiment_id ${experiment_id} \\
        --process_type ${process_type} \\
        --process_id ${process_id} \\
        --srx_accession ${accession} \\
        --organism "${organism}"
    """
}

process PROCESS_TRACKER_FINISH {
    label "tracker_env"
    
    input:
    val accession
    val process_type
    val process_id
    val status

    script:
    def experiment_id = "${accession}"
    """
    python3 ${projectDir}/bin/process_tracker_finish.py \\
        --experiment_id ${experiment_id} \\
        --process_type ${process_type} \\
        --process_id ${process_id} \\
        --status ${status}
    """
}

