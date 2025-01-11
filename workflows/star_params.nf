include { joinReads; } from '../lib/utils.groovy'
include { makeParamSets; validateRequiredColumns; loadBarcodes; loadStarIndices; expandStarParams } from '../lib/star_params.groovy'

// Workflow to run STAR alignment on scRNA-seq data
workflow STAR_PARAMS_WF{
    take:
    ch_accessions
    ch_sra_stat

    main:
    //-- Download subset of reads reads --//    
    // Run prefetch & fastq-dump; only testing on subset of SRR accessions, if many
    ch_fqdump = FASTQ_DUMP(ch_accessions.randomSample(params.max_accessions, 23482))

    // Merge dump logs
    FQDUMP_LOG_MERGE(ch_fqdump.log.collect())
    
    // Join R1 and R2 channels, which will filter out empty R2 records
    ch_fastq = joinReads(ch_fqdump.R1, ch_fqdump.R2)

    // Get read lengths
    SEQKIT_STATS(ch_fastq)

    //-- STAR param search on subsampled reads --//

    // Load barcodes file
    ch_barcodes = loadBarcodes(params)

    // Load star indices
    ch_star_indices = loadStarIndices(params)

    // Pairwise combine samples with barcodes, strand, and star index
    ch_params = makeParamSets(ch_fastq, ch_barcodes, ch_star_indices)

    // Run STAR on subsampled reads, for all pairwise parameter combinations
    STAR_PARAM_SEARCH(ch_params)

    // Format the STAR parameters into a CSV file
    STAR_FORMAT_PARAMS(STAR_PARAM_SEARCH.out)
    
    // Get best parameters
    ch_params_all = STAR_FORMAT_PARAMS.out
        .groupTuple(by: [0,1])
        .join(SEQKIT_STATS.out, by: [0,1])
        .join(ch_sra_stat, by: [0,1])
    STAR_SELECT_PARAMS(ch_params_all)

    // Filter empty params
    ch_star_params_json = STAR_SELECT_PARAMS.out.json
        .filter { sample, accession, json_file -> 
            if(json_file.size() < 5) {
                println "Warning: No valid STAR parameters found for ${sample}; skipping"
            }
            return json_file.size() > 5
        }

    // Extract selected parameters from the JSON files
    ch_star_params = expandStarParams(ch_fastq, ch_star_params_json)

    // Merge STAR parameters to a single set for each sample
    def majorityRule = { list ->
        list.groupBy { it }
            .collectEntries { k, v -> [(k): v.size()] }
            .max { it.value }
            .key
    }

    ch_star_params = ch_star_params.groupTuple()
        .map { sample, accessions, barcodes, star_indices, cell_barcode_lengths, umi_lengths, strands ->
            def combined_params = [barcodes, star_indices, cell_barcode_lengths, umi_lengths, strands].transpose()
            [sample] + majorityRule(combined_params)
        }

    ch_star_params.view{ sample,barcodes_file,star_index,cb_len,umi_len,strand -> 
        def barcodes_name = barcodes_file.tokenize("/").last().tokenize(".")[0]
        def star_index_name = star_index.tokenize("/").last().tokenize(".")[0]
        "Selected parameter set for ${sample}:\n - barcodes: ${barcodes_name}\n - STAR index: ${star_index_name}\n - Barcode length: ${cb_len}\n - UMI length: ${umi_len}\n - Strand: ${strand}\n"
    }

    // Save the final STAR parameters to the database
    STAR_SAVE_FINAL_PARAMS(ch_star_params)

    emit:
    star_params = ch_star_params
}

// Set STAR parameters based on valid barcodes
def saveAsFinalParams(sample, filename) {
    if (filename.endsWith(".csv") || filename.endsWith(".json")){
        filename = filename.tokenize("/").last()
        return "STAR/${sample}/${filename}"
    }
    return null
}

process STAR_SAVE_FINAL_PARAMS {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsFinalParams(sample, filename) }
    label "star_env"

    input:
    tuple val(sample), val(barcodes), val(star_index), val(cell_barcode_length), val(umi_length), val(strand)

    output:
    path "star_params.csv"

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    upload-final-star-params.py \\
      --sample ${sample} \\
      --barcodes ${barcodes} \\
      --star-index ${star_index} \\
      --cell-barcode-length ${cell_barcode_length} \\
      --umi-length ${umi_length} \\
      --strand ${strand}
    """

    stub:
    """
    touch star_params.json
    """
}


process STAR_SELECT_FINAL_PARAMS {
    publishDir file(params.output_dir) / "STAR", mode: "copy", overwrite: true
    label "star_env"

    input:
    tuple val(sample), val(accession), val(metadata),
        path("input*_R1.fastq"), path("input*_R2.fastq"), 
        path(barcodes_file), path(star_index),
        val(cell_barcode_length), val(umi_length), val(strand)

    output:
    path "report.csv"
    
    script:
    """
    create-star-params-report.py star_params_*.json
    """

    stub:
    """
    touch merged_star_params.csv
    """
}

process STAR_MERGE_PARAMS {
    publishDir file(params.output_dir) / "STAR", mode: "copy", overwrite: true
    label "star_env"

    input:
    path "star_params_*.csv"

    output:
    path "merged_star_params.csv"
    
    script:
    """
    csv-merge.py --outfile merged_star_params.csv star_params_*.csv
    """

    stub:
    """
    touch merged_star_params.csv
    """
}

// Set STAR parameters based on valid barcodes
def saveAsParams(sample, accession, filename) {
    if (filename.endsWith(".csv") || filename.endsWith(".json")){
        filename = filename.tokenize("/").last()
        return "STAR/${sample}/${accession}/${filename}"
    }
    return null
}

process STAR_SELECT_PARAMS {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsParams(sample, accession, filename) }
    label "star_env"

    input:
    tuple val(sample), val(accession), path("star_params*.csv"), path(read_stats), path(sra_stats)

    output:
    tuple val(sample), val(accession), path("results/merged_star_params.csv"),    emit: csv
    tuple val(sample), val(accession), path("results/selected_star_params.json"), emit: json
    path "results/select-star-params_log.csv",  emit: "log"
    
    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"
    
    select-star-params.py \\
      --sample ${sample} \\
      --accession ${accession} \\
      $read_stats $sra_stats star_params*.csv
    """

    stub:
    """
    touch star_params.json
    """
}

process STAR_FORMAT_PARAMS {
    label "star_env"

    input:
    tuple val(sample), val(accession), val(metadata), val(params), path(star_summary) 

    output:
    tuple val(sample), val(accession), path("star_params.csv")

    script:
    """
    format-star-params.py \\
      --sample ${params.sample} \\
      --accession ${params.accession} \\
      --strand ${params.strand} \\
      --barcodes-name ${params.barcodes_name} \\
      --barcodes-file ${params.barcodes_file} \\
      --cell-barcode-length ${params.cell_barcode_length} \\
      --umi-length ${params.umi_length} \\
      --organism ${params.organism} \\
      --star-index ${params.star_index} \\
      $star_summary 
    """
}

// Set STAR parameters based on valid barcodes
def saveAsValid(sample, filename) {
    return "param_search/${filename}"
}

// Run STAR alignment on subsampled reads with various parameters to determine which parameters produce the most valid barcodes
process STAR_PARAM_SEARCH {
    label "star_env"
    label "process_medium"

    input:
    tuple val(sample), val(accession), val(metadata), path(fastq_1), path(fastq_2), path(barcodes_file), path(star_index), val(params)

    output:
    tuple val(sample), val(accession), val(metadata), val(params), path("star_summary.csv")

    script:
    """
    # run STAR
    STAR \\
      --readFilesIn $fastq_2 $fastq_1 \\
      --runThreadN ${task.cpus} \\
      --genomeDir ${star_index} \\
      --soloCBwhitelist ${barcodes_file} \\
      --soloCBlen ${params.cell_barcode_length} \\
      --soloUMIlen ${params.umi_length} \\
      --soloStrand ${params.strand} \\
      --soloType CB_UMI_Simple \\
      --clipAdapterType CellRanger4 \\
      --outFilterScoreMin 30 \\
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
      --soloCellFilter EmptyDrops_CR \\
      --soloUMIfiltering MultiGeneUMI_CR \\
      --soloUMIdedup 1MM_CR \\
      --soloFeatures GeneFull \\
      --soloMultiMappers EM \\
      --outSAMtype None \\
      --soloBarcodeReadLength 0 \\
      --outFileNamePrefix results 
    
    # rename output
    mv -f resultsSolo.out/GeneFull/Summary.csv star_summary.csv
    """

    stub:
    """
    touch star_summary.csv
    """
}

// Get read lengths via `seqkit stats`
process SEQKIT_STATS {
    label "download_env"
    label "process_low"

    input:
    tuple val(sample), val(accession), val(metadata), path(fastq_1), path(fastq_2)

    output:
    tuple val(sample), val(accession), path("${sample}_${accession}_stats.tsv")

    script:
    """
    seqkit -j $task.cpus stats -T $fastq_1 $fastq_2 > ${sample}_${accession}_stats.tsv
    """

    stub:
    """
    touch ${sample}_stats.tsv
    """
}

process FQDUMP_LOG_MERGE {
    publishDir file(params.output_dir) / "logs", mode: "copy", overwrite: true
    label "download_env"

    input:
    path "*_log.csv"

    output:
    path "fq-dump.csv"

    script:
    """
    csv-merge.py --outfile fq-dump.csv *_log.csv
    """

    stub:
    """
    touch fq-dump.csv 
    """
}

process FASTQ_DUMP {
    label "download_env"

    input:
    tuple val(sample), val(accession), val(metadata), val(sra_file_size_gb)

    output:
    tuple val(sample), val(accession), val(metadata), path("reads/read_1.fastq"), emit: "R1"
    tuple val(sample), val(accession), val(metadata), path("reads/read_2.fastq"), emit: "R2", optional: true
    path "reads/fq-dump_log.csv", emit: "log"

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    fq-dump.py \\
      --sample ${sample} \\
      --accession ${accession} \\
      --threads ${task.cpus} \\
      --bufsize 10MB \\
      --curcache 50MB \\
      --mem 5GB \\
      --temp TMP_FILES \\
      --min-read-length ${params.min_read_len} \\
      --maxSpotId ${params.max_spots} \\
      --outdir reads \\
      ${accession}

    # remove the temporary files
    rm -rf TMP_FILES
    """

    stub:
    """
    mkdir -p reads
    touch reads/${accession}_1.fastq reads/${accession}_2.fastq
    """
}

process PREFETCH_LOG_MERGE{
    publishDir file(params.output_dir) / "logs", mode: "copy", overwrite: true
    label "download_env"

    input:
    path "*_log.csv"

    output:
    path "prefetch.csv"

    script:
    """
    csv-merge.py --outfile prefetch.csv *_log.csv
    """

    stub:
    """
    touch prefetch.csv 
    """
}
