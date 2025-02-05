workflow { 
    FIND_MTX()

    // list target MTX files
    mtx_files = FIND_MTX.out
        .splitCsv( header: true )
        .map { row -> 
            tuple( row["batch"], row["srx"], row["matrix_path"] ) 
        }
        .groupTuple()

    //mtx_files.view()

    // aggregate mtx files as h5ad
    MTX2H5AD( mtx_files )
}

process MTX2H5AD {
    label "process_high"

    input:
    tuple val(batch), val(srx), val(mtx_path)

    output:
    path "data.h5ad"

    script:
    """
    mtx2h5ad.py --threads ${task.cpus} --srx "$srx" --path "$mtx_path"
    """
}

process FIND_MTX {   
    output:
    path "mtx_files.csv"

    script:
    """
    find-mtx.py \\
      --max-datasets ${params.max_datasets} \\
      --batch-size ${params.batch_size} \\
      --db-uri ${params.db_uri} \\
      ${params.input_dir}
    """
}