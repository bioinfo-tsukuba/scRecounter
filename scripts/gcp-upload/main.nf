workflow { 
    // find target MTX files to add to the database
    FIND_MTX()

    // list target MTX files
    mtx_files = FIND_MTX.out.csv
      .splitCsv( header: true )
      .map{ row -> 
        tuple(row.srx, file(row.matrix_path), file(row.features_path), file(row.barcodes_path))
      }

    // convert to h5ad and publish
    MTX_TO_H5AD( mtx_files, Channel.fromPath(params.tissue_categories) )

    // write parquet after all MTX_TO_H5AD jobs complete
    DB_TO_PARQUET( MTX_TO_H5AD.out.h5ad.collect() )
}

process DB_TO_PARQUET {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, pattern: "metadata/*/*.parquet"
    publishDir file(params.log_dir), mode: "copy", overwrite: true, pattern: "*.log"
    label "process_low"

    input:
    path h5ad_files

    output:
    path "metadata/*/*.parquet", emit: parquet
    path "db-to-parquet.log", emit: log

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    db-to-parquet.py 2>&1 | tee db-to-parquet.log
    """
}

process MTX_TO_H5AD {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, pattern: "h5ad/*/*.h5ad"
    publishDir file(params.log_dir), mode: "copy", overwrite: true, pattern: "*.log"
    label "process_high"
    maxForks 10

    input:
    tuple val(srx), path(mtx_path), path(features_path), path(barcodes_path)
    each path(tissue_categories)

    output:
    path "h5ad/*/${srx}.h5ad",          emit: h5ad
    path "mtx-to-h5ad_${srx}.log", emit: log

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    mtx-to-h5ad.py \\
      --missing-metadata "${params.missing_metadata}" \\
      --tissue-categories "${tissue_categories}" \\
      --srx "${srx}" \\
      --matrix "${mtx_path}" \\
      --publish-path "${params.output_dir}" \\
      2>&1 | tee mtx-to-h5ad_${srx}.log
    """
}

process FIND_MTX {
    publishDir file(params.log_dir), mode: "copy", overwrite: true, pattern: "*.log"
    label "process_low"

    output:
    path "mtx_files.csv", emit: csv
    path "find-mtx.log",  emit: log

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    find-mtx.py \\
      --feature-type ${params.feature_type} \\
      --max-datasets ${params.max_datasets} \\
      ${params.input_dir} \\
      2>&1 | tee find-mtx.log
    """
}