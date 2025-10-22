// Subworkflows
include { DB_ACC_WF } from './workflows/db_acc.nf'
include { STAR_PARAMS_WF } from './workflows/star_params.nf'
include { STAR_FULL_WF } from './workflows/star_full.nf'
include { SRA_STAT } from './lib/utils.nf'
// Process tracker processes
include { PROCESS_TRACKER_START; PROCESS_TRACKER_FINISH } from './workflows/process_tracker.nf'
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

    // run sra-stat on accessions
    ch_sra_stat = SRA_STAT(ch_accessions)
    ch_accessions = addStats(ch_accessions, ch_sra_stat)

    // filter out any accessions with max SRA file size greater than the user-specified size
    ch_accessions = ch_accessions.filter { it[4] <= params.max_sra_size }
    
    // determine best STAR parameters on a subset of reads
    ch_star_params = STAR_PARAMS_WF(ch_accessions, ch_sra_stat)

    // run STAR on all reads with selected parameters
    if (! params.define && ! params.params_only){
        ch_star_results = STAR_FULL_WF(ch_accessions, ch_star_params)
    }

    // Collect final results - all accessions that made it this far are successful
    // Accessions that failed earlier would have been filtered out by Nextflow
    ch_final_results = ch_accessions.map { accession -> [accession[1], "success"] }

    // Single PROCESS_TRACKER_FINISH call with final status
    PROCESS_TRACKER_FINISH(
        ch_final_results.map { accession, status -> accession },
        process_type, 
        process_id, 
        ch_final_results.map { accession, status -> status == "error" ? 1 : 0 }
    )
}

// On complete
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
