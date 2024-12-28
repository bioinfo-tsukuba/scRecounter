include { joinReads; } from '../lib/utils.groovy'

// Workflow to run STAR alignment on scRNA-seq data
workflow STAR_FULL_WF{
    take:
    ch_accessions
    ch_star_params
    
    main:
    //-- Download all reads --//

    // fasterq-dump
    ch_fastq = FASTERQ_DUMP(ch_accessions, ch_star_params.map{ [it[0], it[1]] })
    ch_fastq = joinReads(ch_fastq.R1, ch_fastq.R2)

    // combine reads and star params
    ch_fastq = ch_fastq.join(ch_star_params, by: [0,1])

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
    tuple val(sample), val(accession), path("gene_summary.csv")
    tuple val(sample), val(accession), path("gene_full_summary.csv")
    tuple val(sample), val(accession), path("gene_ex50_summary.csv")
    tuple val(sample), val(accession), path("gene_ex_int_summary.csv")
    tuple val(sample), val(accession), path("velocyto_summary.csv")

    output:
    tuple val(sample), val(accession), path("Summary.csv")

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    star-summary.py \\
      --sample ${sample} \\
      --accession ${accession} \\
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
    tuple val(sample), val(accession), val(metadata),
          path("input*_R1.fastq"), path("input*_R2.fastq"), 
          path(barcodes_file), path(star_index),
          val(cell_barcode_length), val(umi_length), val(strand)

    output:
    tuple val(sample), val(accession), path("resultsSolo.out/Gene/raw/*"),                          emit: gene_raw
    tuple val(sample), val(accession), path("resultsSolo.out/Gene/filtered/*"),                     emit: gene_filt, optional: true
    tuple val(sample), val(accession), path("resultsSolo.out/Gene/Summary.csv"),                    emit: gene_summary
    tuple val(sample), val(accession), path("resultsSolo.out/GeneFull/raw/*"),                      emit: gene_full_raw
    tuple val(sample), val(accession), path("resultsSolo.out/GeneFull/filtered/*"),                 emit: gene_full_filt, optional: true
    tuple val(sample), val(accession), path("resultsSolo.out/GeneFull/Summary.csv"),                emit: gene_full_summary
    tuple val(sample), val(accession), path("resultsSolo.out/GeneFull_Ex50pAS/raw/*"),              emit: gene_ex50_raw
    tuple val(sample), val(accession), path("resultsSolo.out/GeneFull_Ex50pAS/filtered/*"),         emit: gene_ex50_filt, optional: true
    tuple val(sample), val(accession), path("resultsSolo.out/GeneFull_Ex50pAS/Summary.csv"),        emit: gene_ex50_summary
    tuple val(sample), val(accession), path("resultsSolo.out/GeneFull_ExonOverIntron/raw/*"),       emit: gene_ex_int_raw
    tuple val(sample), val(accession), path("resultsSolo.out/GeneFull_ExonOverIntron/filtered/*"),  emit: gene_ex_int_filt, optional: true
    tuple val(sample), val(accession), path("resultsSolo.out/GeneFull_ExonOverIntron/Summary.csv"), emit: gene_ex_int_summary
    tuple val(sample), val(accession), path("resultsSolo.out/Velocyto/raw/*"),                      emit: velocyto_raw
    tuple val(sample), val(accession), path("resultsSolo.out/Velocyto/filtered/*"),                 emit: velocyto_filt, optional: true
    tuple val(sample), val(accession), path("resultsSolo.out/Velocyto/Summary.csv"),                emit: velocyto_summary

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
    """
}

def saveAsSTAR(sample, accession, filename) {
    if (filename.endsWith(".mtx") || filename.endsWith(".tsv") || filename.endsWith(".csv")) {
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
    memory { (16.GB + Math.round(sra_file_size_gb).GB) * task.attempt }
    time { (4.h + (sra_file_size_gb * 1.1).h) * task.attempt }
    disk 750.GB, type: "local-ssd"
    cpus 8
    maxRetries 2

    input:
    tuple val(sample), val(accession), val(metadata), val(sra_file_size_gb)
    tuple val(sample), val(accession)

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