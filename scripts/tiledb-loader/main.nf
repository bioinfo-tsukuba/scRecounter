workflow { 
    // find target MTX files to add to the database
    FIND_MTX()

    // list target MTX files
    mtx_files = FIND_MTX.out.csv
        .splitCsv( header: true )
        .map { row -> 
            tuple( row["batch"], row["srx"], row["matrix_path"] ) 
        }
        .groupTuple()

    // aggregate mtx files as h5ad
    MTX_TO_H5AD( mtx_files )

    // add the h5ad files to the database
    H5AD_TO_DB( MTX_TO_H5AD.out.h5ad.buffer( size: params.h5ad_batch_size, remainder: true ) )
}

process H5AD_TO_DB {
    publishDir file(params.log_dir), mode: "copy", overwrite: true
    label "process_high"
    maxForks 1

    input:
    path "?.h5ad"

    output:
    path "h5ad_to_db.log", emit: log

    script:
    """
    h5ad-to-db.py \\
      --threads ${task.cpus} \\
      --db-uri ${params.db_uri} \\
      *.h5ad 2>&1 | tee h5ad_to_db.log
    """
}

process MTX_TO_H5AD {
    publishDir file(params.log_dir) , mode: "copy", overwrite: true, pattern: "*.log"
    label "process_high"
    maxForks 6

    input:
    tuple val(batch), val(srx), val(mtx_path)

    output:
    path "data.h5ad",                      emit: h5ad
    path "mtx_to_h5ad_batch-${batch}.log", emit: log

    script:
    """
    mtx-to-h5ad.py \\
      --threads ${task.cpus} \\
      --missing-metadata "${params.missing_metadata}" \\
      --srx "$srx" \\
      --path "$mtx_path" \\
      2>&1 | tee mtx_to_h5ad_batch-${batch}.log
    """
}

process FIND_MTX {
    publishDir file(params.log_dir), mode: "copy", overwrite: true, pattern: "*.log"
    label "process_low"

    output:
    path "mtx_files.csv", emit: csv
    path "find_mtx.log",  emit: log

    script:
    """
    find-mtx.py \\
      --feature-type ${params.feature_type} \\
      --max-datasets ${params.max_datasets} \\
      --batch-size ${params.mtx_batch_size} \\
      --db-uri ${params.db_uri} \\
      ${params.input_dir} \\
      2>&1 | tee find_mtx.log
    """
}