import groovy.json.JsonSlurper

def expandStarParams(ch_fastq, ch_star_params_json) {
    return ch_fastq.join(ch_star_params_json, by: [0,1])
        .map{ sample, accession, read1, read2, json_file -> 
            def params = new JsonSlurper().parseText(json_file.text)
            def barcodes_file = params.barcodes_file
            def star_index = params.star_index
            def cell_barcode_length = params.cell_barcode_length
            def umi_length = params.umi_length
            def strand = params.strand
            return [sample, accession, 
                    read1, read2, barcodes_file, star_index, 
                    cell_barcode_length, umi_length, strand]
        }
}


def makeParamSets(ch_subsample, ch_barcodes, ch_star_indices) {
    return ch_subsample
        .combine(Channel.of("Forward", "Reverse"))
        .combine(ch_barcodes)
        .combine(ch_star_indices)
        .map { sample, accession, r1, r2, strand, barcodes_name, cb_len, umi_len, barcodes_file, organism, star_index ->
            def params = [
                sample: sample,
                accession: accession,
                strand: strand,
                barcodes_name: barcodes_name,
                cell_barcode_length: cb_len,
                umi_length: umi_len,
                barcodes_file: barcodes_file,
                organism: organism,
                star_index: star_index
            ]
            return [sample, accession, r1, r2, barcodes_file, star_index, params] 
        }
}


def validateRequiredColumns(row, required) {
    def missing = required.findAll { !row.containsKey(it) }
    if (missing) {
        error "Missing columns in the input CSV file: ${missing}"
    }
}

def loadBarcodes(params) {
    return Channel
        .fromPath(params.barcodes, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def req_columns = ["name", "cell_barcode_length", "umi_length", "file_path"]
            validateRequiredColumns(row, req_columns)
            // remove special characters
            row.name = row.name.replaceAll("\\s", "_")
            return [
                row.name, 
                row.cell_barcode_length.toInteger(), 
                row.umi_length.toInteger(),
                row.file_path
            ]
        }
}

def loadStarIndices(params) {
    return Channel
        .fromPath(params.star_indices, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def req_columns = ["organism", "star_index"]
            validateRequiredColumns(row, req_columns)
            // remove special characters
            row.organism = row.organism.replaceAll("\\s", "_")
            return [row.organism, row.star_index]
        }
}