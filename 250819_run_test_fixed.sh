#!/bin/bash
outdir="../results_test_URL_fixed"

if [ ! -d "$outdir" ]; then
  mkdir -p "$outdir"
fi

# Run the Nextflow pipeline with the fixed CSV format
nextflow -log $outdir/run_test.log \
    run main.nf  \
    -work-dir /data01/kariyayama/CellIO/tmp/work_by_URL_fixed/    \
    -c config/resources.config \
    -profile conda,trace,report,vm,vm_dev,dev   \
    --barcodes data/my_barcodes.csv     \
    --star_indices data/my_star_indices.csv     \
    --accessions data/250819_accession_by_URL_fixed.csv     \
    --output_dir $outdir # \
#     --reads_with_barcodes_cutoff 0.01 \
#     --max_spots 2000000
