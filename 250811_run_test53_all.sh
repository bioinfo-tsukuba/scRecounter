# Run the Nextflow pipeline with the specified parameters
for target_acc in 250811_test53_001 250811_test53_002 250811_test53_003 250811_test53_004 250811_test53_005;
do
    # 出力ファイル
    outdir="${HOME}/CellIO/250811_test53/${target_acc}"    
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
        --accessions   ${HOME}/CellIO/Test53_target_list/${target_acc}.csv \
        --output_dir   $outdir
    
    # .fastqを削除
    ./remove_recursive.sh $cur_tmp fastq
    # .sraを削除
    ./remove_recursive.sh $cur_tmp sra
done

