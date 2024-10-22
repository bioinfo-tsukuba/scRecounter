// Subworkflows
include { DOWNLOAD_WF } from './workflows/download.nf'
include { READS_WF } from './workflows/reads.nf'
include { READ_QC_WF } from './workflows/read_qc.nf'
include { STAR_WF } from './workflows/star.nf'

// Main workflow
workflow {
    if (params.accessions){
        // DOWNLOAD
        DOWNLOAD_WF()
    } else {
        // Load reads 
        READS_WF()
    }

    // READ_QC
    //READ_QC_WF(READS_WF.out.fastq)

    // STAR
    //STAR_WF(READS_WF.out.fastq)
}

// On complete
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
