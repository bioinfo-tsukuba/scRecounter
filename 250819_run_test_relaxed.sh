#\!/bin/bash
outdir="../results_test_URL"

if [ \! -d "$outdir" ]; then
  mkdir -p "$outdir"
fi

# Run the Nextflow pipeline with relaxed threshold
nextflow -log $outdir/run_test_relaxed.log \
    run main.nf  \
    -work-dir /data01/kariyayama/CellIO/tmp/work_by_URL_4/    \
    -c config/resources.config \
    -profile conda,trace,report,vm,vm_dev,dev   \
    --barcodes data/my_barcodes.csv     \
    --star_indices data/my_star_indices.csv     \
    --accessions data/250819_accession_by_URL_2.csv     \
    --output_dir $outdir \
    --reads_with_barcodes_cutoff 0.02 \
    --max_spots 2000000
