workflow DOWNLOAD_WF {
    main:
    // Load accessions from file
    ch_accessions = Channel
        .fromPath(params.accessions, checkIfExists: true)
        .splitCsv(header: false, sep: ",")
        .map { row -> row[0] }

    // Get vdb-dump info
    //VDB_DUMP_INFO(ch_accessions)
    //VDB_DUMP_MERGE(VDB_DUMP_INFO.out.collect{ it[1] })

    // prefetch
    PREFETCH(ch_accessions)
    
    // fasterq-dump
    FASTERQ_DUMP(PREFETCH.out)

    // Merge by sample


}

process FASTERQ_DUMP {
    conda "envs/download.yml"
    label "process_low"
    //scratch true

    input:
    tuple val(accession), path(sra_file)

    output:
    tuple val(accession), path("reads/${accession}_1.fastq"), emit: "R1"
    tuple val(accession), path("reads/${accession}_2.fastq"), emit: "R2", optional: true

    script:
    """
    # Run fasterq-dump
    fasterq-dump \\
      --threads ${task.cpus} \\
      --bufsize 5MB \\
      --curcache 50MB \\
      --mem 5GB \\
      --temp TMP_FILES \\
      --outdir reads \\
      --split-files \\
      --force \\
      ${sra_file}
    
    # Remove the temporary files and sra file
    #rm -rf TMP_FILES ${sra_file}
    """

    stub:
    """
    mkdir -p reads
    touch reads/${accession}_1.fastq reads/${accession}_2.fastq
    """
}

process PREFETCH {
    conda "envs/download.yml"
    label "process_low"

    input:
    val accession

    output:
    tuple val(accession), path("prefetch_out/${accession}/${accession}.sra")

    script:
    """
    prefetch.py \\
      --max-size 5000 \\
      --outdir prefetch_out \\
      ${accession}
    """

    stub:
    """
    mkdir -p prefetch_out/${accession}
    touch prefetch_out/${accession}/${accession}.sra
    """
}

process VDB_DUMP_MERGE {
    publishDir file(params.outdir) / "download", mode: "copy", overwrite: true
    conda "envs/download.yml"
    label "process_low"

    input:
    path info

    output:
    path "info.csv"

    script:
    """
    vdb-dump_merge.py $info > info.csv
    """
}

process VDB_DUMP_INFO {
    conda "envs/download.yml"
    label "process_low"

    input:
    val accession

    output:
    tuple val(accession), path("${accession}_info.txt")

    script:
    """
    vdb-dump --info ${accession} > ${accession}_info.txt
    """
}


// --list-sizes
// --ascp-path
// vdb-validate

/*
vdb-config --cfg 
vdb-config --set "/repository/user/main/public/root=/scratch/multiomics/nickyoungblut/sra/"
vdb-config --set "/repository/user/cache/root=/scratch/multiomics/nickyoungblut/sra_cache/"
prefetch --transport http --max-size 1G --output-directory prefetch_out SRR13711615
prefetch --max-size 1G --output-directory prefetch_out SRR13711615



fasterq-dump \\
      --threads ${task.cpus} \\
      --bufsize 10MB \\
      --curcache 50MB
      --mem 5GB \\
      --temp TMP_FILES \\
      --outdir data \\
      --split-files \\
      --spots 10000 \\
      --force \\
      ${accession} 
    
    # remove the temporary files
    rm -rf TMP_FILES
*/

/*
process FASTQ_DUMP {
    conda "envs/download.yml"
    label "process_low"

    input:
    val accession

    output:
    tuple val(accession), path("data/${accession}_1.fastq"), emit: "R1"
    tuple val(accession), path("data/${accession}_2.fastq"), emit: "R2", optional: true

    script:
    """
    # Run fastq-dump
    fastq-dump \\
      --outdir data \\
      --split-files \\
      --maxSpotId 10000 \\
      ${accession} 
    """
}
*/