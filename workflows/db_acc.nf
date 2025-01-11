include { readAccessions } from '../lib/download.groovy'

workflow DB_ACC_WF {
    main:
    // obtain accessions from the database
    ch_accessions = GET_DB_ACCESSIONS()
    ch_accessions.ifEmpty { println 'No accessions found in the scRecounter database' }

    emit:
    ch_accessions
}


process GET_DB_ACCESSIONS {
    label "download_env"

    output:
    path "accessions.csv"

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    get-db-accessions.py --max-srx ${params.max_samples}
    """
}