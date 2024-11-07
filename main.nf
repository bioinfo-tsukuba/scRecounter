// Subworkflows
include { DOWNLOAD_WF } from './workflows/download.nf'
include { READS_WF } from './workflows/reads.nf'
include { READ_QC_WF } from './workflows/read_qc.nf'
include { STAR_PARAMS_WF } from './workflows/star_params.nf'
include { STAR_FULL_WF } from './workflows/star_full.nf'

// Main workflow
workflow {
    if (params.star_params){
        // User-provided selected STAR parameters
        ch_fastq = Channel
            .fromPath(params.star_params, checkIfExists: true)
            .splitCsv(header: true, sep: ',')
            .map { row ->
                def req_columns = ["sample", "fastq_1", "fastq_2", "barcodes_file", "star_index", 
                                   "cell_barcode_length", "umi_length", "strand"]
                def miss_columns = req_columns.findAll { !row.containsKey(it) }
                if (miss_columns) {
                error "Missing columns in the input CSV file: ${miss_columns}"
            }
            // remove special characters from the sample name
            row.sample = row.sample.replaceAll("\\s", "_")
            return [row.sample, row.fastq_1, row.fastq_2, row.barcodes_file, row.star_index, 
                    row.cell_barcode_length, row.umi_length, row.strand]
        }
    } else {
        // Obtain the STAR parameters 
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

        // Select STAR parameters
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
