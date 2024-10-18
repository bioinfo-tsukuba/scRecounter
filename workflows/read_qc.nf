workflow READ_QC_WF{
  take:
  ch_fastq

  main:
  // seqkit stats
  SEQKIT_STATS(ch_fastq)
  SEQKIT_STATS_MERGE(SEQKIT_STATS.out.collect())
  
}

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

process SEQKIT_STATS {
  conda "envs/read_qc.yml"
  label "process_low"

  input:
  tuple val(sample), path(fastq_1), path(fastq_2)

  output:
  path "${sample}.tsv"

  script:
  """
  seqkit -j $task.cpus stats -a -T $fastq_1 $fastq_2 > ${sample}.tsv
  """

  stub:
  """
  touch ${sample}.tsv
  """
}