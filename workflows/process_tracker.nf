// Process tracker processes
process PROCESS_TRACKER_START {
    label "tracker_env"
    
    input:
    val accession
    val process_type
    val process_id

    script:
    def experiment_id = "${accession}_${workflow.start.format('yyyyMMdd_HHmmss')}"
    """
    python3 ${projectDir}/bin/process_tracker_start.py \\
        --experiment_id ${experiment_id} \\
        --process_type ${process_type} \\
        --process_id ${process_id}
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
    def experiment_id = "${accession}_${workflow.start.format('yyyyMMdd_HHmmss')}"
    """
    python3 ${projectDir}/bin/process_tracker_finish.py \\
        --experiment_id ${experiment_id} \\
        --process_type ${process_type} \\
        --process_id ${process_id} \\
        --status ${status}
    """
}

