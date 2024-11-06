// Subworkflows
include { DOWNLOAD_WF } from './workflows/download.nf'
include { READS_WF } from './workflows/reads.nf'
include { READ_QC_WF } from './workflows/read_qc.nf'
include { STAR_PARAMS_WF } from './workflows/star_params.nf'
include { STAR_FULL_WF } from './workflows/star_full.nf'

// Main workflow
workflow {
    if (params.accessions){
        // Download reads
        DOWNLOAD_WF()
        ch_fastq = DOWNLOAD_WF.out.fastq
        ch_sra_stat = DOWNLOAD_WF.out.sra_stat
    } else {
        // Load existing reads
        ch_fastq = READS_WF()
        ch_sra_stat = Channel.empty() // TODO: update
    }

    // QC on reads
    //READ_QC_WF(ch_fastq)

    // STAR
    STAR_PARAMS_WF(ch_fastq, ch_sra_stat)
    STAR_FULL_WF(STAR_PARAMS_WF.out)
}

// On complete
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
