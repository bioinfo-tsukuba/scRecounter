// Subworkflows
include { DB_ACC_WF } from './workflows/db_acc.nf'
include { DOWNLOAD_WF } from './workflows/download.nf'
include { READS_WF } from './workflows/reads.nf'
include { READ_QC_WF } from './workflows/read_qc.nf'
include { STAR_PARAMS_WF } from './workflows/star_params.nf'
include { STAR_FULL_WF } from './workflows/star_full.nf'
// utils
include { readStarParams } from './lib/utils.groovy'

// Main workflow
workflow {    
    if (params.star_params){
        // User-provided selected STAR parameters
        println("Using provided STAR parameters.")
        ch_fastq = readStarParams(params.star_params)
    } else {
        if (params.fastq){
            // Load existing reads
            println("Using provided fastq files.")
            ch_fastq = READS_WF()
            ch_sra_stat = Channel.empty() // TODO: update
        } else {
            if (params.accessions == "" || params.accessions == true) {
                // Obtain accessions from SRA
                println("No accessions provided. Accessions will be obtained from SRA.")
                ch_accessions = DB_ACC_WF()
            } else {
                // Use the provided accessions
                println("Using provided accessions.")
                ch_accessions = Channel.fromPath(params.accessions, checkIfExists: true)
            }
            // Download reads from SRA
            DOWNLOAD_WF(ch_accessions)
            ch_fastq = DOWNLOAD_WF.out.fastq
            ch_sra_stat = DOWNLOAD_WF.out.sra_stat
        }
        // determine STAR parameters
        ch_fastq = STAR_PARAMS_WF(ch_fastq, ch_sra_stat)
    } 

    // Run STAR on all reads with selected parameters
    if (! params.define){
        STAR_FULL_WF(ch_fastq)
    }
}

// On complete
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
