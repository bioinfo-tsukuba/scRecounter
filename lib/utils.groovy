def readStarParams(star_params_file){
    return Channel
        .fromPath(star_params_file, checkIfExists: true)
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
}
