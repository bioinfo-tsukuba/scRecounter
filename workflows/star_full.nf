include { joinReads; } from '../lib/utils.groovy'

// Workflow to run STAR alignment on scRNA-seq data
workflow STAR_FULL_WF{
    take:
    ch_accessions
    ch_star_params
    
    main:
    //-- Download all reads --//
    // filter out samples that lack a set of selected parameters
    samples_to_keep = ch_star_params.map{ it[0] }.collect()
    ch_accessions_filt = ch_accessions.filter{
        samples_to_keep.val.contains(it[0])
    }

    // fasterq-dump to download all reads
    ch_fastq = FASTERQ_DUMP(ch_accessions_filt)
    ch_fastq = joinReads(ch_fastq.R1, ch_fastq.R2)

    // combine reads and star params
    ch_fastq = ch_fastq.map{ sample, accession, metadata, R1, R2 -> [sample, R1, R2] }
        .groupTuple()
        .join(ch_star_params)

    //-- Run STAR with the selected parameters on all reads --//
    // run STAR
    STAR_FULL(ch_fastq)

    // summarize the STAR results
    STAR_FULL_SUMMARY(
        STAR_FULL.out.gene_summary,
        STAR_FULL.out.gene_full_summary,
        STAR_FULL.out.gene_ex50_summary,
        STAR_FULL.out.gene_ex_int_summary,
        STAR_FULL.out.velocyto_summary
    )
}

process STAR_FULL_SUMMARY {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsSTAR(sample, accession, filename) }
    label "star_env"
    label "process_high"

    input:
    tuple val(sample), path("gene_summary.csv")
    tuple val(sample), path("gene_full_summary.csv")
    tuple val(sample), path("gene_ex50_summary.csv")
    tuple val(sample), path("gene_ex_int_summary.csv")
    tuple val(sample), path("velocyto_summary.csv")

    output:
    tuple val(sample), path("Summary.csv")

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    star-summary.py \\
      --sample ${sample} \\
      gene_summary.csv \\
      gene_full_summary.csv \\
      gene_ex50_summary.csv \\
      gene_ex_int_summary.csv \\
      velocyto_summary.csv
    """
}

process STAR_FULL {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsSTAR(sample, accession, filename) }
    label "star_env"
    label "process_high"

    input:
    tuple val(sample), path("input*_R1.fastq"), path("input*_R2.fastq"), 
          path(barcodes_file), path(star_index),
          val(cell_barcode_length), val(umi_length), val(strand)

    output: 
    tuple val(sample), path("resultsSolo.out/Gene/Summary.csv"),                    emit: gene_summary
    tuple val(sample), path("resultsSolo.out/GeneFull/Summary.csv"),                emit: gene_full_summary
    tuple val(sample), path("resultsSolo.out/GeneFull_Ex50pAS/Summary.csv"),        emit: gene_ex50_summary
    tuple val(sample), path("resultsSolo.out/GeneFull_ExonOverIntron/Summary.csv"), emit: gene_ex_int_summary
    tuple val(sample), path("resultsSolo.out/Velocyto/Summary.csv"),                emit: velocyto_summary
    tuple val(sample), path("resultsSolo.out/*/raw/*"),                             emit: raw
    tuple val(sample), path("resultsSolo.out/*/filtered/*"),                        emit: filt, optional: true
    tuple val(sample), path("resultsSolo.out/*/*.stats.gz"),                        emit: stats, optional: true
    tuple val(sample), path("resultsSolo.out/*/*.txt.gz"),                          emit: txt, optional: true

    script:
    """
    # run STAR
    R1=\$(printf "%s," input*_R1.fastq)
    R1=\${R1%,} 
    R2=\$(printf "%s," input*_R2.fastq)
    R2=\${R2%,}
    STAR \\
      --readFilesIn \$R2 \$R1 \\
      --runThreadN ${task.cpus} \\
      --genomeDir ${star_index} \\
      --soloCBwhitelist ${barcodes_file} \\
      --soloUMIlen ${umi_length} \\
      --soloStrand ${strand} \\
      --soloCBlen ${cell_barcode_length} \\
      --soloType CB_UMI_Simple \\
      --clipAdapterType CellRanger4 \\
      --outFilterScoreMin 30 \\
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
      --soloCellFilter EmptyDrops_CR \\
      --soloUMIfiltering MultiGeneUMI_CR \\
      --soloUMIdedup 1MM_CR \\
      --soloFeatures Gene GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS Velocyto \\
      --soloMultiMappers EM Uniform \\
      --outSAMtype None \\
      --soloBarcodeReadLength 0 \\
      --outFileNamePrefix results

    # gzip the results
    mkdir -p resultsSolo.out
    find resultsSolo.out -type f -name "*.stats" | xargs -P ${task.cpus} gzip
    find resultsSolo.out -type f -name "*.txt" | xargs -P ${task.cpus} gzip
    find resultsSolo.out -type f -name "*.tsv" | xargs -P ${task.cpus} gzip
    find resultsSolo.out -type f -name "*.mtx" | xargs -P ${task.cpus} gzip
    """
}

def saveAsSTAR(sample, accession, filename) {
    def extensions = [".mtx.gz", ".tsv.gz", ".txt.gz", ".stats.gz", ".csv"]
    if (extensions.any { filename.endsWith(it) }) {
        def parts = filename.tokenize("/")
        if (parts.size() > 1) {
            return "STAR/${sample}/${accession}/" + parts[1..-1].join('/')
        } else {
            return "STAR/${sample}/${accession}/" + parts[0]
        }
    } 
    return null
}

process FASTERQ_DUMP {
    label "download_env"
    memory { (16.GB + Math.round(sra_file_size_gb / 2).GB) * task.attempt }
    time { (4.h + (sra_file_size_gb * 1.1).h) * task.attempt }
    disk 750.GB, type: "local-ssd"
    cpus 8
    maxRetries 2

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
      --bufsize 200MB \\
      --curcache 1GB \\
      --mem 12GB \\
      --temp TMP_FILES \\
      --max-size-gb ${params.max_sra_size} \\
      --min-read-length ${params.min_read_len} \\
      --maxSpotId ${params.max_spots} \\
      --outdir reads \\
      ${accession}
    """

    stub:
    """
    mkdir -p reads
    touch reads/${accession}_1.fastq reads/${accession}_2.fastq
    """
}