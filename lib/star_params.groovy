import groovy.json.JsonSlurper

def expandStarParams(ch_fastq, ch_star_params_json) {
    ch_params = ch_fastq.join(ch_star_params_json, by: [0,1])
        .map{ sample, accession, metadata, read1, read2, json_file -> 
            def params = new JsonSlurper().parseText(json_file.text)
            def barcodes_file = params.barcodes_file
            def star_index = params.star_index
            def cell_barcode_length = params.cell_barcode_length
            def umi_length = params.umi_length
            def strand = params.strand
            return [sample, accession,
                    barcodes_file, star_index, 
                    cell_barcode_length, umi_length, strand]
        }

    // status on number of parameter combinations
    ch_params
        .ifEmpty("No valid parameter set selected")
        .view{ sample,accession,barcodes_file,star_index,cb_len,umi_len,strand -> 
            def barcodes_name = barcodes_file.tokenize("/").last().tokenize(".")[0]
            def star_index_name = star_index.tokenize("/").last().tokenize(".")[0]
            "Selected parameter set:\n - barcodes: ${barcodes_name}\n - STAR index: ${star_index_name}\n - Barcode length: ${cb_len}\n - UMI length: ${umi_len}\n - Strand: ${strand}\n"
        }
    return ch_params
}

def makeParamSets(ch_subsample, ch_barcodes, ch_star_indices) {
    ch_params = ch_subsample
        .combine(Channel.of("Forward", "Reverse"))
        .combine(ch_barcodes)
        .combine(ch_star_indices)
        .map { sample, accession, metadata, r1, r2, strand, barcodes_name, cb_len, umi_len, barcodes_file, organism, star_index ->
            if (metadata["organism"] != organism) {
                return null
            }
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
            return [sample, accession, metadata, r1, r2, barcodes_file, star_index, params] 
        }
        .filter { it != null }

    // status on number of parameter combinations
    ch_params
        .ifEmpty("No valid parameter set found")
        .count().view{ count -> "Number of parameter sets to test: ${count}" }
    return ch_params
}

def validateRequiredColumns(row, required) {
    def missing = required.findAll { !row.containsKey(it) }
    if (missing) {
        error "Missing columns in the input CSV file: ${missing}"
    }
}

def loadBarcodes(params) {
    ch_barcodes = Channel
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
    // status on number of barcodes
    ch_barcodes
        .ifEmpty("No barcodes found in the input CSV file")
        .count().view{ count -> "Number of input barcodes: ${count}" }
    return ch_barcodes
}

def loadStarIndices(params) {
    ch_indices = Channel
        .fromPath(params.star_indices, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def req_columns = ["organism", "star_index"]
            validateRequiredColumns(row, req_columns)
            // remove special characters
            row.organism = row.organism.replaceAll("\\s", "_")
            return [row.organism, row.star_index]
        }
    // status on number of star indices
    ch_indices
        .ifEmpty("No star indices found in the input CSV file")
        .count().view{ count -> "Number of input star indices: ${count}" }
    return ch_indices
}