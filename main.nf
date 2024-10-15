// Subworkflows
include { DOWNLOAD_WF } from './workflows/download.nf'

// Main workflow
workflow {

}

// On complete
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
