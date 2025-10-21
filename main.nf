// Subworkflows
include { DB_ACC_WF } from './workflows/db_acc.nf'
include { STAR_PARAMS_WF } from './workflows/star_params.nf'
include { STAR_FULL_WF } from './workflows/star_full.nf'
include { SRA_STAT } from './lib/utils.nf'
// util functions
include { readAccessions; addStats; } from './lib/utils.groovy'

// Main workflow
workflow { 
    // PROCESS TRACKER INTEGRATION - COMMENTED OUT
    // Initialize ProcessTracker for experiment tracking
    // experiment_id = "${workflow.runName}_${workflow.start.format('yyyyMMdd_HHmmss')}"
    // process_type = "scRecounter"
    // process_id = "version_0.1"
    // 
    // // Start process tracking
    // """
    // python3 ${projectDir}/bin/process_tracker_start.py \\
    //     --experiment_id ${experiment_id} \\
    //     --process_type ${process_type} \\
    //     --process_id ${process_id}
    // """
    
    if (params.accessions == "" || params.accessions == true) {
        // Obtain accessions from SRA
        println "No accessions provided. Accessions will be obtained from SRA."
        ch_accessions = DB_ACC_WF()
    } else {
        // Use the provided accessions
        println "Using provided accessions."
        ch_accessions = Channel.fromPath(params.accessions, checkIfExists: true)
    }

    // read accessions file
    ch_accessions = readAccessions(ch_accessions)

    // run sra-stat on accessions
    ch_sra_stat = SRA_STAT(ch_accessions)
    ch_accessions = addStats(ch_accessions, ch_sra_stat)

    // filter out any accessions with max SRA file size greater than the user-specified size
    ch_accessions = ch_accessions.filter { it[4] <= params.max_sra_size }
    
    // determine best STAR parameters on a subset of reads
    ch_star_params = STAR_PARAMS_WF(ch_accessions, ch_sra_stat)

    // run STAR on all reads with selected parameters
    if (! params.define){
        STAR_FULL_WF(ch_accessions, ch_star_params)
    }
}

// On complete
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    
    // PROCESS TRACKER INTEGRATION - COMMENTED OUT
    // Finish process tracking with comprehensive error reporting
    // def status = workflow.success ? 0 : 1
    // def error_message = ""
    // if (!workflow.success) {
    //     // Collect error information from multiple sources
    //     def errorSources = []
    //     if (workflow.errorMessage) errorSources.add("Workflow: ${workflow.errorMessage}")
    //     if (workflow.errorReport) errorSources.add("Report: ${workflow.errorReport}")
    //     error_message = errorSources.join(" | ")
    // }
    // 
    // """
    // python3 ${projectDir}/bin/process_tracker_finish.py \\
    //     --experiment_id ${experiment_id} \\
    //     --process_type ${process_type} \\
    //     --process_id ${process_id} \\
    //     --status ${status} \\
    //     --error_message "${error_message}"
    // """
}
