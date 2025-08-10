# outdir="results_test"
outdir="results_test53_head20"


if [ ! -d "$outdir" ]; then
  mkdir -p "$outdir"
fi

# Run the Nextflow pipeline with the specified parameters
nextflow -log $outdir/run_test.log \
    run main.nf  \
    -work-dir /data01/kariyayama/CellIO/tmp/work/    \
    -c config/resources.config \
    -profile conda,trace,report,vm,vm_dev,dev   \
    --barcodes data/my_barcodes.csv     \
    --star_indices data/my_star_indices.csv     \
    --accessions data/SraRunInfo_test53_selected_head20.csv     \
    --output_dir $outdir \
    --resume
      
