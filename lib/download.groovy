def readAccessions(accessions_input){
    // Read the input CSV file with the sample names and SRA accessions
    ch_acc = accessions_input
        .splitCsv(header: true, sep: ",")
        .map { row ->
            def req_columns = ["sample", "accession"]
            def miss_columns = req_columns.findAll { !row.containsKey(it) }
            if (miss_columns) {
                error "Missing columns in the input CSV file: ${miss_columns}"
            }
            // remove special characters from the sample name
            row.sample = row.sample.replaceAll("\\s", "_")
            def result = [row.sample, row.accession]
            // add optional, metadata columns
            def metadata = [:]
            ["organism", "tech_10x"].each { col ->
                metadata[col] = row.containsKey(col) ? row[col] : ""
            }
            result << metadata
            return result
        }

    // print srx values
    ch_acc
        .map{ sample, accession, metadata -> sample }
        .distinct()
        .collect() 
        .map{ it.join(',') }
        .view{ "SRX accessions: ${it}" } 

    return ch_acc
}

/*
def addStats(ch_accessions, ch_sra_stat){
    // Add sra-stats channel to the accessions channel
    ch_stats = ch_sra_stat
        .map{ sample,acc,csv -> csv }
        .splitCsv(header: true, sep: ",")
        .map{ row -> [row.accession, row.file_size_gb.toDouble()] }
        .join(ch_sra_stat.map{ sample,acc,csv -> [acc, sample] }, by: [0])
        .map{ acc, size, sample -> [sample, acc, size] }
    return ch_accessions.join(ch_stats, by: [0,1]) // sample, acc, metadata, size
}

def joinReads(ch_read1, ch_read2){
    // Join the read1 and read2 channels
    ch_metadata = ch_read1.map{ sample,accession,metadata,fastq -> [sample,accession,metadata] }

    return ch_read1
        .map{ sample,accession,metadata,fastq -> [sample,accession,fastq] }
        .join(
            ch_read2.map{ sample,accession,metadata,fastq -> [sample,accession,fastq] }, 
            by: [0,1]
        )
        .join(
            ch_metadata, by: [0,1]
        )
        .map{ 
            sample,accession,fastq1,fastq2,metadata -> [sample,accession,metadata,fastq1,fastq2] 
        }
}
*/
