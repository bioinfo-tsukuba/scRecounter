include { saveAsLog } from '../lib/utils.groovy'
include { readAccessions } from '../lib/download.groovy'

workflow DB_ACC_WF {
    main:
    // obtain accessions from the database
    ch_accessions = GET_DB_ACCESSIONS()
    ch_accessions.csv.ifEmpty { println 'No accessions found in the scRecounter database' }

    emit:
    ch_accessions.csv
}


process GET_DB_ACCESSIONS {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename) }
    label "download_env"

    output:
    path "accessions.csv",      emit: "csv"
    path "${task.process}.log", emit: "log"

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    get-db-accessions.py \\
      --max-srx ${params.max_samples} \\
      2>&1 | tee ${task.process}.log
    """
}

