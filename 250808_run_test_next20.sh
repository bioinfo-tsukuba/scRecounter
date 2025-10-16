# outdir="results_test"
outdir="results_test53_next20"

if [ ! -d "$outdir" ]; then
  mkdir -p "$outdir"
fi

# Run the Nextflow pipeline with the specified parameters
nextflow -log $outdir/run_test.log \
    run main.nf  \
    -work-dir /data01/kariyayama/CellIO/tmp/work_next20    \
    -c config/resources.config \
    -profile conda,trace,report,vm,vm_dev,dev   \
    --barcodes data/my_barcodes.csv     \
    --star_indices data/my_star_indices.csv     \
    --accessions data/SraRunInfo_test53_selected_next20.csv     \
    --output_dir $outdir

