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
workflow{
    if (params.accessions == "" || params.accessions == true) {
        // Obtain accessions from SRA
        println "No accessions provided. Accessions will be obtained from SRA."
        // ch_accessions = DB_ACC_WF()
    } else {
        // Use the provided accessions
        println "Using provided accessions."
        ch_accessions = Channel
            .fromPath(params.accessions, checkIfExists: true)
            .tap { dbg_fromPath }
        dbg_fromPath.view { v -> println "Before readAccessions ${v}" }
    }

    // read accessions file
    ch_accessions = readAccessions(ch_accessions)
        .tap { dbg_readAccessions }
    dbg_readAccessions.view { v -> println "After readAccessions ${v}"  }

    // run sra-stat on accessions
    ch_sra_stat = SRA_STAT(ch_accessions)
        .tap { dbg_sra_stat }
    dbg_sra_stat.view { v -> println "After SRA_STAT ${v}"  }

    // add sra-stat results to accessions channel
    ch_accessions = addStats(ch_accessions, ch_sra_stat)
        .tap { dbg_addStats }
    dbg_addStats.view { v -> println "After addStats ${v}" }
}

