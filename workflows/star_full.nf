include { joinReads; saveAsLog; } from '../lib/utils.groovy'

// Workflow to run STAR alignment on scRNA-seq data
workflow STAR_FULL_WF{
    take:
    ch_accessions
    ch_star_params
    
    main:
    //-- Download all reads --//
    // filter out samples that lack a set of selected parameters
    ch_accessions_filt = ch_accessions.combine(
        ch_star_params.map{ it[0] }.unique(), by: 0
    )

    // fasterq-dump to download all reads
    ch_fastq = FASTERQ_DUMP(ch_accessions_filt)
    ch_fastq = joinReads(ch_fastq.R1, ch_fastq.R2)

    // For accessions lacking paired reads from fasterq-dump, fallback to fastq-dump
    ch_accessions_fallback = ch_accessions_filt
        .join(
            ch_fastq.map{ it -> [it[0], it[1], true] }, 
            by: [0,1],
            remainder: true
        )
        .filter{ it -> it[5] != true }
        .map{ it -> it[0..4] }

    // run fastq-dump on the fallback accessions
    ch_fastq_fallback = FASTQ_DUMP(ch_accessions_fallback)
    ch_fastq_fallback = joinReads(ch_fastq_fallback.R1, ch_fastq_fallback.R2)
    ch_fastq_fallback.count().view{ count -> "No. of fastq-dump fallback accessions: $count" }

    // combine the fasterq-dump and fastq-dump results
    ch_fastq = ch_fastq.mix(ch_fastq_fallback)
    ch_fastq.count().view{ count -> "No. of fast(er)q-dump accessions: $count" }

    // 個別accessionデータを保持（groupTuple前）
    ch_individual_accessions = ch_fastq
        .map{ sample, accession, metadata, R1, R2 -> [sample, accession, metadata] }

    // combine reads and star params
    ch_fastq = ch_fastq
        .map{ sample, accession, metadata, R1, R2 -> [sample, R1, R2] }
        .groupTuple()
        .join(ch_star_params)

    //-- Run STAR with the selected parameters on all reads --//
    // run STAR
    STAR_FULL(ch_fastq)

    emit:
    // 個別accessionの結果（成功・失敗両方を追跡）
    individual_results = STAR_FULL.out.status
        .join(ch_individual_accessions.groupTuple(by: 0))  // sampleでグループ化されたaccessionとjoin
        .flatMap { sample, exit_status, accessions, metadatas ->
            // exit_statusに基づいて全accessionのステータスを設定
            def results = []
            def status = (exit_status as Integer) == 0 ? 0 : 1
            for (int i = 0; i < accessions.size(); i++) {
                results << [sample, accessions[i], status]  // [sample, accession, status]
            }
            return results
        }
    
    // 既存互換性のためのsample集約結果（成功のみ）
    success_results = individual_results
        .filter { sample, accession, status -> status == 0 }
        .groupTuple(by: 0)
        .map { sample, accessions, statuses -> 
            [sample, accessions[0], 0]  // 代表accessionを使用
        }
}

process STAR_FULL {
    tag "${sample}_${accession}"
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsSTAR(sample, filename) }
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample) }
    label "star_env"
    label "process_high"
    errorStrategy 'ignore'  // エラーでもワークフローを継続
    disk { [request: (375 * (task.attempt > 1 ? 2 : 1)).GB, type: 'local-ssd'] }
    machineType { 
        def options = ['n2-*', 'n2d-*']
        return options[new Random().nextInt(options.size())]
    }

    input:
    tuple val(sample), path("input*_R1.fastq.gz"), path("input*_R2.fastq.gz"), 
          path(barcodes_file), path(star_index),
          val(cell_barcode_length), val(umi_length), val(strand)

    output: 
    tuple val(sample), path("resultsSolo.out/Gene/Summary.csv"),                    emit: gene_summary, optional: true
    tuple val(sample), path("resultsSolo.out/GeneFull/Summary.csv"),                emit: gene_full_summary, optional: true
    tuple val(sample), path("resultsSolo.out/GeneFull_Ex50pAS/Summary.csv"),        emit: gene_ex50_summary, optional: true
    tuple val(sample), path("resultsSolo.out/GeneFull_ExonOverIntron/Summary.csv"), emit: gene_ex_int_summary, optional: true
    tuple val(sample), path("resultsSolo.out/Velocyto/Summary.csv"),                emit: velocyto_summary, optional: true
    tuple val(sample), path("resultsSolo.out/*/raw/*"),                             emit: raw, optional: true
    tuple val(sample), path("resultsSolo.out/*/filtered/*"),                        emit: filt, optional: true
    tuple val(sample), path("resultsSolo.out/*/*.stats.gz"),                        emit: stats, optional: true
    tuple val(sample), path("resultsSolo.out/*/*.txt.gz"),                          emit: txt, optional: true
    tuple val(sample), env(EXIT_STATUS),                                            emit: status
    path "${task.process}.log",                                                     emit: "log"

    script:
    """
    echo "Running STAR for ${sample}" > ${task.process}.log
    
    # Initialize EXIT_STATUS
    EXIT_STATUS=0

    R1=\$(printf "%s," input*_R1.fastq.gz)
    R1=\${R1%,} 
    R2=\$(printf "%s," input*_R2.fastq.gz)
    R2=\${R2%,}
    STAR \\
      --readFilesIn \$R2 \$R1 \\
      --readFilesCommand zcat \\
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
      --outFileNamePrefix results \\
      2>&1 | tee -a ${task.process}.log
    
    # Capture STAR exit status
    EXIT_STATUS=\$?
    
    if [ \$EXIT_STATUS -eq 0 ]; then
        # gzip the results only on success
        mkdir -p resultsSolo.out
        find resultsSolo.out -type f -name "*.stats" | xargs -P ${task.cpus} gzip
        find resultsSolo.out -type f -name "*.txt" | xargs -P ${task.cpus} gzip
        find resultsSolo.out -type f -name "*.tsv" | xargs -P ${task.cpus} gzip
        find resultsSolo.out -type f -name "*.mtx" | xargs -P ${task.cpus} gzip
    else
        # Create empty output directories on failure
        mkdir -p resultsSolo.out/Gene resultsSolo.out/GeneFull resultsSolo.out/GeneFull_Ex50pAS resultsSolo.out/GeneFull_ExonOverIntron resultsSolo.out/Velocyto
        echo "STAR failed with exit code \$EXIT_STATUS" >> ${task.process}.log
    fi
    """
}

def saveAsSTAR(sample, filename) {
    def extensions = [".mtx.gz", ".tsv.gz", ".txt.gz", ".stats.gz", ".csv"]
    if (extensions.any { filename.endsWith(it) }) {
        def parts = filename.tokenize("/")
        if (parts.size() > 1) {
            return "STAR/${sample}/" + parts[1..-1].join('/')
        } else {
            return "STAR/${sample}/" + parts[0]
        }
    } 
    return null
}

process FASTQ_DUMP {
    tag "${sample}_${accession}"
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample, accession) }
    label "download_env"
    maxRetries 1
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    cpus 4
    memory { 4.GB * task.attempt }
    time { (6.h + (sra_file_size_gb * 0.8).h) * task.attempt }
    disk {[request: 375.GB, type: 'local-ssd']}
    machineType { 
        def options = ['n2-*', 'c2-*', 'n2d-*', 'c2d-*']
        return options[new Random().nextInt(options.size())]
    }
    
    input:
    tuple val(sample), val(accession), val(download_url), val(metadata), val(sra_file_size_gb)

    output:
    tuple val(sample), val(accession), val(metadata), path("reads/read_1.fastq.gz"), emit: "R1"
    tuple val(sample), val(accession), val(metadata), path("reads/read_2.fastq.gz"), emit: "R2", optional: true
    path "${task.process}.log",                                                   emit: "log"

    script:
    def sra_input = download_url ?: accession
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    echo "Downloading ${accession} for ${sample}" > ${task.process}.log
    echo "sra-stat file size: ${sra_file_size_gb} GB" >> ${task.process}.log

    echo "Running fastq-dump as backup for fasterq-dump" >> ${task.process}.log
    fq-dump.py \\
      --sample ${sample} \\
      --accession ${accession} \\
      --threads ${task.cpus} \\
      --min-read-length ${params.min_read_len} \\
      --outdir reads \\
      --maxSpotId ${params.fallback_max_spots} \\
      ${sra_input} \\
      2>&1 | tee -a ${task.process}.log
    """

    stub:
    """
    mkdir -p reads
    touch reads/read_1.fastq.gz reads/read_2.fastq.gz ${task.process}.log
    """
}

process FASTERQ_DUMP {
    tag "${sample}_${accession}"
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample, accession) }
    label "download_env"
    maxRetries 1
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    cpus 6
    memory { 16.GB * task.attempt }
    time { (10.h + (sra_file_size_gb * 0.8).h) * task.attempt }
    disk { 
        def disk_size = 
            sra_file_size_gb > 260 ? 2625.GB :
            sra_file_size_gb > 220 ? 2250.GB :
            sra_file_size_gb > 170 ? 1875.GB :
            sra_file_size_gb > 120 ? 1500.GB :
            sra_file_size_gb > 60 ? 1125.GB :
            sra_file_size_gb > 30 ? 750.GB :
            375.GB
        disk_size = disk_size + (375 * (task.attempt - 1)).GB
        [request: disk_size, type: 'local-ssd'] 
    }
    machineType { 
        def options = ['n2-*', 'c2-*', 'n2d-*', 'c2d-*'] //, 'n1-*']
        return options[new Random().nextInt(options.size())]
    }
    
    input:
    tuple val(sample), val(accession), val(download_url), val(metadata), val(sra_file_size_gb)

    output:
    tuple val(sample), val(accession), val(metadata), path("reads/read_1.fastq.gz"), emit: "R1", optional: true
    tuple val(sample), val(accession), val(metadata), path("reads/read_2.fastq.gz"), emit: "R2", optional: true
    path "${task.process}.log",                                                   emit: "log"

    script:
    def sra_input = download_url ?: accession
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    echo "Downloading ${accession} for ${sample}" > ${task.process}.log
    echo "sra-stat file size: ${sra_file_size_gb} GB" >> ${task.process}.log

    # run prefetch and fasterq-dump
    fq-dump.py \\
      --sample ${sample} \\
      --accession ${accession} \\
      --threads ${task.cpus} \\
      --bufsize 200MB \\
      --curcache 1GB \\
      --mem 12GB \\
      --max-size-gb ${params.max_sra_size} \\
      --min-read-length ${params.min_read_len} \\
      --temp ${params.fasterq_tmp} \\
      --outdir reads \\
      ${sra_input} \\
      2>&1 | tee -a ${task.process}.log
    """

    stub:
    """
    mkdir -p reads
    touch reads/read_1.fastq.gz reads/read_2.fastq.gz ${task.process}.log
    """
}
