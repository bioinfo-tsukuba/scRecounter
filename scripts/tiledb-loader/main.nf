workflow { 
    // find target MTX files to add to the database
    FIND_MTX()

    // list target MTX files
    mtx_files = FIND_MTX.out
        .splitCsv( header: true )
        .map { row -> 
            tuple( row["batch"], row["srx"], row["matrix_path"] ) 
        }
        .groupTuple()

    // aggregate mtx files as h5ad
    MTX_TO_H5AD( mtx_files )

    // add the h5ad files to the database
    H5AD_TO_DB( MTX_TO_H5AD.out.buffer( size: params.h5ad_batch_size, remainder: true ) )
}

process H5AD_TO_DB {
    label "process_high"

    input:
    path "?.h5ad"

    script:
    """
    h5ad-to-db.py --db-uri ${params.db_uri} *.h5ad
    """
}

process MTX_TO_H5AD {
    label "process_high"

    input:
    tuple val(batch), val(srx), val(mtx_path)

    output:
    path "data.h5ad"

    script:
    """
    mtx-to-h5ad.py --threads ${task.cpus} --srx "$srx" --path "$mtx_path"
    """
}

process FIND_MTX {   
    output:
    path "mtx_files.csv"

    script:
    """
    find-mtx.py \\
      --max-datasets ${params.max_datasets} \\
      --batch-size ${params.mtx_batch_size} \\
      --db-uri ${params.db_uri} \\
      ${params.input_dir}
    """
}