// Subworkflows
//include { DOWNLOAD_WF } from './workflows/download.nf'
include { INPUT_WF } from './workflows/input.nf'

// Main workflow
workflow {
    // load input
    INPUT_WF()
    // 
    
}

// On complete
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
