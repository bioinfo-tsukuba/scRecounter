// Process tracker processes
process PROCESS_TRACKER_START {
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

// Error-aware wrapper processes
process SRA_STAT_WITH_TRACKING {
    errorStrategy 'ignore'
    
    input:
    val accession
    val process_type
    val process_id

    output:
    tuple val(accession), val("success"), emit: success, optional: true
    tuple val(accession), val("error"), emit: error, optional: true

    script:
    """
    set +e
    # Run SRA_STAT process logic here or call existing process
    # For now, this is a placeholder - you'll need to integrate with actual SRA_STAT
    echo "Running SRA_STAT for ${accession}"
    
    # Check if process succeeded
    if [ \$? -eq 0 ]; then
        echo "success" > status.txt
    else
        echo "error" > status.txt
    fi
    """
}

process STAR_PARAMS_WF_WITH_TRACKING {
    errorStrategy 'ignore'
    
    input:
    val accession
    val process_type
    val process_id

    output:
    tuple val(accession), val("success"), emit: success, optional: true
    tuple val(accession), val("error"), emit: error, optional: true

    script:
    """
    set +e
    echo "Running STAR_PARAMS_WF for ${accession}"
    
    # Check if process succeeded
    if [ \$? -eq 0 ]; then
        echo "success" > status.txt
    else
        echo "error" > status.txt
    fi
    """
}

process STAR_FULL_WF_WITH_TRACKING {
    errorStrategy 'ignore'
    
    input:
    val accession
    val process_type
    val process_id

    output:
    tuple val(accession), val("success"), emit: success, optional: true
    tuple val(accession), val("error"), emit: error, optional: true

    script:
    """
    set +e
    echo "Running STAR_FULL_WF for ${accession}"
    
    # Check if process succeeded
    if [ \$? -eq 0 ]; then
        echo "success" > status.txt
    else
        echo "error" > status.txt
    fi
    """
}