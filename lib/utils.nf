process SRA_STAT {
    tag "${sample}_${accession}"
    label "download_env"
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    disk 10.GB

    input:
    tuple val(sample), val(accession), val(download_url), val(metadata)

    output:
    tuple val(sample), val(accession), path("sra-stat.csv")

    afterScript """
    # Read the main script exit status from Nextflow's .exitcode file
    if [ -f .exitcode ]; then
        EXIT_CODE=\$(cat .exitcode)
    else
        EXIT_CODE=0
    fi
    echo "=== AFTERSCRIPT DEBUG START ==="
    echo "AfterScript EXIT_CODE: \$EXIT_CODE"
    echo "Sample: ${sample}"
    echo "Accession: ${accession}" 
    echo "Experiment ID: ${sample}_${accession}"
    echo "Project Dir: ${projectDir}"
    echo "Checking if python script exists..."
    ls -la ${projectDir}/bin/process_tracker_finish.py || echo "Script not found!"
    
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "=== CALLING PROCESS_TRACKER_FINISH ==="
        echo "SRA_STAT failed with exit code \$EXIT_CODE for ${sample}_${accession}"
        python3 ${projectDir}/bin/process_tracker_finish.py \\
            --experiment_id "${sample}_${accession}" \\
            --process_type "scRecounter" \\
            --process_id "version_0.1" \\
            --status 2 \\
            --error_message "SRA_STAT process ignored due to error (exit code: \$EXIT_CODE)" || echo "PROCESS_TRACKER_FINISH failed"
    else
        echo "=== SRA_STAT SUCCESS, NO TRACKER CALL NEEDED ==="
    fi
    echo "=== AFTERSCRIPT DEBUG END ==="
    """

    script:
    """
    echo "=== SRA_STAT DEBUG START ==="
    echo "Sample: ${sample}"
    echo "Accession: ${accession}"
    echo "Download URL: ${download_url}"
    echo "Metadata: ${metadata}"
    echo "Project Dir: ${projectDir}"
    echo "Task Process: ${task.process}"
    echo "Task Work Dir: ${task.workDir}"
    echo "=== EXECUTING sra-stat.py ==="
    
    sra-stat.py ${accession}
    SCRIPT_EXIT_CODE=\$?
    
    echo "=== SRA_STAT SCRIPT COMPLETED ==="
    echo "Script exit code: \$SCRIPT_EXIT_CODE"
    echo "=== SRA_STAT DEBUG END ==="
    
    exit \$SCRIPT_EXIT_CODE
    """

    stub:
    """
    touch sra-stat.csv
    """
}