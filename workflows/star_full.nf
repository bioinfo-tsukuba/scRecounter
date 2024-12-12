// Workflow to run STAR alignment on scRNA-seq data
workflow STAR_FULL_WF{
    take:
    ch_fastq

    main:
    //-- Run STAR with the selected parameters --//
    STAR_FULL(ch_fastq)
}

// STAR alignment with all reads and selected parameters
def saveAsSTAR(sample, accession, filename) {
    if (filename.endsWith(".mtx") || filename.endsWith(".tsv") || filename.endsWith(".csv")) {
        def parts = filename.tokenize("/")
        return "STAR/${sample}/${accession}/" + parts[1..-1].join('/')
    } 
    return null
}

process STAR_FULL {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsSTAR(sample, accession, filename) }
    container "us-east1-docker.pkg.dev/c-tc-429521/sc-recounter-star/sc-recounter-star:0.1.0"
    conda "envs/star.yml"
    label "process_high"

    input:
    tuple val(sample), val(accession), 
          path("input*_R1.fastq"), path("input*_R2.fastq"), 
          path(barcodes_file), path(star_index),
          val(cell_barcode_length), val(umi_length), val(strand)

    output:
    tuple val(sample), path("resultsSolo.out/Gene/raw/*"),                          emit: gene_raw
    tuple val(sample), path("resultsSolo.out/Gene/filtered/*"),                     emit: gene_filt, optional: true
    tuple val(sample), path("resultsSolo.out/Gene/Summary.csv"),                    emit: gene_summary
    tuple val(sample), path("resultsSolo.out/GeneFull/raw/*"),                      emit: gene_full_raw
    tuple val(sample), path("resultsSolo.out/GeneFull/filtered/*"),                 emit: gene_full_filt, optional: true
    tuple val(sample), path("resultsSolo.out/GeneFull/Summary.csv"),                emit: gene_full_summary
    tuple val(sample), path("resultsSolo.out/GeneFull_Ex50pAS/raw/*"),              emit: gene_ex50_raw
    tuple val(sample), path("resultsSolo.out/GeneFull_Ex50pAS/filtered/*"),         emit: gene_ex50_filt, optional: true
    tuple val(sample), path("resultsSolo.out/GeneFull_Ex50pAS/Summary.csv"),        emit: gene_ex50_summary
    tuple val(sample), path("resultsSolo.out/GeneFull_ExonOverIntron/raw/*"),       emit: gene_ex_int_raw
    tuple val(sample), path("resultsSolo.out/GeneFull_ExonOverIntron/filtered/*"),  emit: gene_ex_int_filt, optional: true
    tuple val(sample), path("resultsSolo.out/GeneFull_ExonOverIntron/Summary.csv"), emit: gene_ex_int_summary
    tuple val(sample), path("resultsSolo.out/Velocyto/raw/*"),                      emit: velocyto_raw
    tuple val(sample), path("resultsSolo.out/Velocyto/filtered/*"),                 emit: velocyto_filt, optional: true
    tuple val(sample), path("resultsSolo.out/Velocyto/Summary.csv"),                emit: velocyto_summary

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
