// Read QC workflow for fastq files
workflow READ_QC_WF{
    take:
    ch_fastq

    main:
    // Flatten the channel to process each read file separately
    ch_fastq_flat = ch_fastq.flatMap { sample, fastq_1, fastq_2 -> 
        [ [sample, 'R1', fastq_1], [sample, 'R2', fastq_2] ] 
    }

    // Run seqkit stats
    SEQKIT_STATS(ch_fastq_flat)
        .collectFile(
          name: "seqkit-stats.tsv", 
          storeDir:  file(params.outdir) / "read_qc", 
          newLine: false, keepHeader: true
        )
    //SEQKIT_STATS_MERGE(SEQKIT_STATS.out.collect())
}

// Run `seqkit stats` on fastq files
process SEQKIT_STATS {
    conda "envs/read_qc.yml"
    label "process_low"

    input:
    tuple val(sample), val(read), path(fastq)

    output:
    path "${sample}_${read}.tsv"

    script:
    """
    seqkit -j $task.cpus stats -a -T $fastq > ${sample}_${read}.tsv
    """

    stub:
    """
    touch ${sample}_${read}.tsv
    """
}

/*
process SEQKIT_STATS_MERGE {
    publishDir file(params.outdir) / "read_qc" , mode: "copy", overwrite: true
    conda "envs/read_qc.yml"

    input:
    path tsv

    output:
    path "seqkit-stats.tsv"

    script:
    """
    seqkit-stats_merge.py $tsv > seqkit-stats.tsv
    """

    stub:
    """
    touch seqkit-stats.tsv
    """
}
*/