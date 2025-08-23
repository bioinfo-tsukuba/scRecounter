# Run the Nextflow pipeline with the specified parameters
target_acc="250823_test_set"
# target_acc="250823_test_set_ddbj"

# 出力ファイル
outdir="${HOME}/CellIO/${target_acc}"    
if [ ! -d "$outdir" ]; then
  mkdir -p "$outdir"
fi

cur_tmp="/data01/kariyayama/CellIO/tmp/work_${target_acc}"
nextflow -log $outdir/${target_acc}.log \
    run main.nf  \
    -work-dir $cur_tmp       \
    -c             config/resources.config  \
    -profile       conda,trace,report,vm,vm_dev,dev     \
    --barcodes     data/my_barcodes.csv     \
    --star_indices data/my_star_indices.csv \
    --accessions   ${target_acc}.csv \
    --output_dir   $outdir

# .fastqを削除
./remove_recursive.sh $cur_tmp fastq.gz
# .sraを削除
./remove_recursive.sh $cur_tmp sra


