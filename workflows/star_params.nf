include { makeParamSets; validateRequiredColumns; loadBarcodes; loadStarIndices; expandStarParams } from '../lib/star_params.groovy'

// Workflow to run STAR alignment on scRNA-seq data
workflow STAR_PARAMS_WF{
    take:
    ch_fastq
    ch_sra_stat

    main:
    //-- for each sample (accession), run STAR on subset of reads with various parameters to determine which params produce the most valid barcodes --/
    
    // Subsample reads
    SUBSAMPLE_READS(ch_fastq)

    // Get read lengths
    SEQKIT_STATS(SUBSAMPLE_READS.out)

    // Load barcodes file
    ch_barcodes = loadBarcodes(params)

    // Load star indices
    ch_star_indices = loadStarIndices(params)

    // Pairwise combine samples with barcodes, strand, and star index
    ch_params = makeParamSets(SUBSAMPLE_READS.out, ch_barcodes, ch_star_indices)

    // Run STAR on subsampled reads, for all pairwise parameter combinations
    STAR_PARAM_SEARCH(ch_params)
    
    // Format the STAR parameters
    STAR_FORMAT_PARAMS(STAR_PARAM_SEARCH.out)

    // Get best parameters
    ch_params_all = STAR_FORMAT_PARAMS.out
        .groupTuple()
        .join(SEQKIT_STATS.out, by: 0)
        .join(ch_sra_stat, by: 0)
    STAR_SELECT_PARAMS(ch_params_all)

    // Merge param CSVs
    STAR_MERGE_PARAMS(STAR_SELECT_PARAMS.out.csv.collect { it[1] })

    // Filter empty params
    ch_star_params_json = STAR_SELECT_PARAMS.out.json
        .filter { sample, json_file -> 
            if(json_file.size() < 5) {
                println "Warning: No valid STAR parameters found for ${sample}; skipping"
            }
            return json_file.size() > 5
        }

    // Extract selected parameters from the JSON files
    ch_fastq = expandStarParams(ch_fastq, ch_star_params_json)

    emit:
    fastq = ch_fastq
}

process STAR_SELECT_PARAMS_REPORT {
    publishDir file(params.outdir) / "STAR", mode: "copy", overwrite: true
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-star/sc-recounter-star:0.1.0"
    conda "envs/star.yml"

    input:
    path "star_params_*.json"

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
    publishDir file(params.outdir) / "STAR", mode: "copy", overwrite: true
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-star/sc-recounter-star:0.1.0"
    conda "envs/star.yml"

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
def saveAsParams(sample, filename) {
    if (filename.endsWith(".csv") || filename.endsWith(".json")){
        filename = filename.tokenize("/").last()
        return "${sample}/${filename}"
    }
    return null
}

process STAR_SELECT_PARAMS {
    publishDir file(params.outdir) / "STAR", mode: "copy", overwrite: true, saveAs: { filename -> saveAsParams(sample, filename) }
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-star/sc-recounter-star:0.1.0"
    conda "envs/star.yml"

    input:
    tuple val(sample), path("star_params*.csv"), path(read_stats), path(sra_stats)

    output:
    tuple val(sample), path("results/merged_star_params.csv"),    emit: csv
    tuple val(sample), path("results/selected_star_params.json"), emit: json
    path "results/select-star-params_log.csv",                    emit: "log"
    
    script:
    """
    select-star-params.py $read_stats $sra_stats star_params*.csv
    """

    stub:
    """
    touch star_params.json
    """
}

process STAR_FORMAT_PARAMS {
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-star/sc-recounter-star:0.1.0"
    conda "envs/star.yml"

    input:
    tuple val(sample), val(params), path(star_summary) 

    output:
    tuple val(sample), path("star_params.csv")

    script:
    """
    format-star-params.py \\
      --sample ${params.sample} \\
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
    return "${sample}/param_search/${filename}"
}

// Run STAR alignment on subsampled reads with various parameters to determine which parameters produce the most valid barcodes
process STAR_PARAM_SEARCH {
    //publishDir file(params.outdir) / "STAR", mode: "copy", overwrite: true, saveAs: { filename -> saveAsValid(sample, filename) }
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-star/sc-recounter-star:0.1.0"
    conda "envs/star.yml"
    label "process_medium"

    input:
    tuple val(sample), path(fastq_1), path(fastq_2), path(barcodes_file), path(star_index), val(params)

    output:
    tuple val(sample), val(params), path("star_summary.csv")

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
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-download/sc-recounter-download:0.1.0"
    conda "envs/read_qc.yml"
    label "process_low"

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)

    output:
    tuple val(sample), path("${sample}_stats.tsv")

    script:
    """
    seqkit -j $task.cpus stats -T $fastq_1 $fastq_2 > ${sample}_stats.tsv
    """

    stub:
    """
    touch ${sample}_stats.tsv
    """
}

// Subsample reads
process SUBSAMPLE_READS {
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-download/sc-recounter-download:0.1.0"
    conda "envs/read_qc.yml"
    label "process_low"

    input:
    tuple val(sample), path("input*_R1.fq"), path("input*_R2.fq")

    output:
    tuple val(sample), path("${sample}_R1.fq"), path("${sample}_R2.fq")

    script: 
    """
    subsample.py --num-seqs ${params.subsample} --out-file ${sample}_R1.fq input*_R1.fq
    subsample.py --num-seqs ${params.subsample} --out-file ${sample}_R2.fq input*_R2.fq
    """
    
    stub:
    """
    touch ${sample}_R1.fq ${sample}_R2.fq
    """
}


