// Subworkflows
include { DB_ACC_WF } from './workflows/db_acc.nf'
include { STAR_PARAMS_WF } from './workflows/star_params.nf'
include { STAR_FULL_WF } from './workflows/star_full.nf'
include { SRA_STAT } from './lib/utils.nf'
// Process tracker processes
include { PROCESS_TRACKER_START; PROCESS_TRACKER_FINISH; SRA_STAT_WITH_TRACKING; STAR_PARAMS_WF_WITH_TRACKING; STAR_FULL_WF_WITH_TRACKING } from './processes/process_tracker.nf'
// util functions
include { readAccessions; addStats; } from './lib/utils.groovy'

// Main workflow
workflow { 
    // Initialize ProcessTracker for experiment tracking
    process_type = "scRecounter"
    process_id = "version_0.1"
    
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

    // Start process tracking for each accession
    PROCESS_TRACKER_START(ch_accessions, process_type, process_id)

    // Run SRA_STAT with error tracking
    ch_sra_results = SRA_STAT_WITH_TRACKING(ch_accessions, process_type, process_id)
    
    // Collect all results for final reporting
    ch_all_results = ch_sra_results.error.mix(ch_sra_results.success)

    // Run STAR_PARAMS_WF with error tracking
    ch_star_params_results = STAR_PARAMS_WF_WITH_TRACKING(
        ch_sra_results.success.map { accession, status -> accession },
        process_type, 
        process_id
    )
    
    // Add STAR_PARAMS results to collection
    ch_all_results = ch_all_results.mix(ch_star_params_results.error).mix(ch_star_params_results.success)

    // Run STAR_FULL_WF with error tracking
    ch_star_full_results = STAR_FULL_WF_WITH_TRACKING(
        ch_star_params_results.success.map { accession, status -> accession },
        process_type, 
        process_id
    )
    
    // Add STAR_FULL results to collection
    ch_all_results = ch_all_results.mix(ch_star_full_results.error).mix(ch_star_full_results.success)

    // Single PROCESS_TRACKER_FINISH call with all results
    PROCESS_TRACKER_FINISH(
        ch_all_results.map { accession, status -> accession },
        process_type, 
        process_id, 
        ch_all_results.map { accession, status -> status == "error" ? 1 : 0 }
    )
}

// On complete
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
