# Run the Nextflow pipeline with the specified parameters
target_file="250825_test_set"

# 出力ファイル
outdir="${HOME}/CellIO/${target_file}"
if [ ! -d "$outdir" ]; then
  mkdir -p "$outdir"
fi

cur_tmp="/data01/kariyayama/CellIO/tmp/work_${target_file}"
nextflow -log $outdir/${target_file}.log \
    run main.nf  \
    -work-dir $cur_tmp       \
    -with-trace "$outdir/${target_file}_trace.txt" \
    -c             config/resources.config  \
    -profile       conda,trace,report,vm,vm_dev,dev     \
    --barcodes     data/my_barcodes.csv     \
    --star_indices data/my_star_indices.csv \
    --accessions   ${target_file}.csv \
    --output_dir   $outdir

# .fastqを削除
./remove_recursive.sh $cur_tmp fastq.gz
# .sraを削除
./remove_recursive.sh $cur_tmp sra

# 結果をTSVで出力
~/miniforge3/envs/cell-io-mappingctl/bin/python \
    ../cell-io-mappingctl/check_screcounter_trace.py \
    --input-csv  $outdir/${target_file}_trace.txt \
    --output-csv $outdir/${target_file}_status.tsv

