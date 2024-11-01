// Subworkflows
include { DOWNLOAD_WF } from './workflows/download.nf'
include { READS_WF } from './workflows/reads.nf'
include { READ_QC_WF } from './workflows/read_qc.nf'
include { STAR_WF } from './workflows/star.nf'

// Main workflow
workflow {
    if (params.accessions){
        // Download reads
        ch_fastq = DOWNLOAD_WF()
    } else {
        // Load existing reads
        ch_fastq = READS_WF()
    }

    // QC on reads
    //READ_QC_WF(ch_fastq)

    // STAR
    STAR_WF(ch_fastq)
}

// On complete
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
