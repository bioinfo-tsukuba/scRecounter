def readAccessions(accessions_input){
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
        .collect() 
        .map{ it.join(',') }
        .view{ "SRX accessions: ${it}" } 

    return ch_acc
}

def joinReads(ch_read1, ch_read2){
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


/*
def joinReads(ch_read1, ch_read2){
    return ch_read1.join(ch_read2, by: [0, 1], remainder: true)
        .filter { sample, accession, r1, r2 -> 
            if(r2 == null) {
                println "Warning: Read 2 is empty for ${sample}; skipping"
            }
            return r2 != null
        }
        .groupTuple()
        .map { sample, accession, fastq_1, fastq_2 ->
            return [sample, fastq_1.flatten(), fastq_2.flatten()]
        }
}
*/