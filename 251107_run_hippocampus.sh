# Run the Nextflow pipeline with the specified parameters
dataset_dir="${HOME}/CellIO/cell-io-mappingctl/data/sra_metadata_split_url"
tmp_dir="/data01/kariyayama/CellIO/tmp"

version=`cat VERSION`
target_label="25_hippocampus_20251104_182024"

# create analysis target list
ls ${dataset_dir}/${target_label}_* | \
    xargs wc -l | \
    awk '{if($1 > 1) print $2}' | \
    grep PAIRED | \
    grep -v .not_found.csv \
    > 251107_analysis_target_${target_label}.txt

# target_acc="26_immune_20251104_182046_human_PAIRED_1028_urls"
cat 251107_analysis_target_${target_label}.txt | \
    awk -F "/" '{print $8}' | \
    sed 's/.csv//g' | \
    while read target_acc; do
        acc_file="${dataset_dir}/${target_acc}.csv"
        outprefix="/cell-io/${version}/${target_acc}"
        cur_tmp="${tmp_dir}/work_${target_acc}"
        
        # 出力ファイル
        if [ ! -d "$outprefix" ]; then
          mkdir -p "$outprefix"
        fi
        
        nextflow -log "${outprefix}.log" \
            run main.nf  \
            -work-dir      $cur_tmp       \
            -c             config/resources.config  \
            -profile       conda,trace,report,vm,vm_dev,dev     \
            --barcodes     data/my_barcodes.csv     \
            --star_indices data/my_star_indices.csv \
            --accessions   ${acc_file} \
            --output_dir   $outprefix
        
        # .fastqを削除
        ./remove_recursive.sh $cur_tmp fastq.gz
        # .sraを削除
        ./remove_recursive.sh $cur_tmp sra
    done
