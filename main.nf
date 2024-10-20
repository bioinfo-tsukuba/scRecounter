// Subworkflows
//include { DOWNLOAD_WF } from './workflows/download.nf'
include { INPUT_WF } from './workflows/input.nf'
include { READ_QC_WF } from './workflows/read_qc.nf'
include { STAR_WF } from './workflows/star.nf'

// Main workflow
workflow {
    // load input
    INPUT_WF()

    // READ_QC
    READ_QC_WF(INPUT_WF.out.fastq)

    // STAR
    STAR_WF(INPUT_WF.out.fastq)
}

// On complete
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
