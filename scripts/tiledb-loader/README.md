tiledb loader
=============

# Production

```bash
cd /home/nickyoungblut/dev/nextflow/scRecounter/scripts/tiledb-loader
```

### scRecounter: prod3

```bash
nextflow run main.nf -profile conda,slurm,report,trace \
  -N nick.youngblut@arcinstitute.org \
  -ansi-log false \
  --max_datasets 1000 \
  --log_dir /scratch/multiomics/nickyoungblut/tiledb-loader/logs/prod3/GeneFull_Ex50pAS \
  --db_uri /scratch/multiomics/nickyoungblut/tiledb-loader/tiledb-soma_GeneFull_Ex50pAS \
  --input_dir /processed_datasets/scRecount/scRecounter/prod3 \
  > logs/prod3-1.log 2>&1
```

Remove working directory

```bash
rm -rf /scratch/multiomics/nickyoungblut/nextflow-work/tiledb-loader/
```

### Chris' Cell-X-Gene data

```bash
nextflow run main.nf -profile conda,slurm,report,trace \
  -N nick.youngblut@arcinstitute.org \
  -ansi-log false \
  --max_datasets 10000 \
  --missing_metadata allow \
  --log_dir /scratch/multiomics/nickyoungblut/tiledb-loader/logs/cellxgene/GeneFull_Ex50pAS \
  --db_uri  /scratch/multiomics/nickyoungblut/tiledb-loader/tiledb-soma_GeneFull_Ex50pAS \
  --input_dir /processed_datasets/scRecount/cellxgene/counted_SRXs \
  > logs/cellxgene-1.log 2>&1
```


# Dev

Local run

```bash
nextflow run main.nf -profile conda,vm,dev -resume
```

Slurm run

```bash
nextflow run main.nf -profile conda,slurm,dev -resume
```

## Test prod

```bash
nextflow run main.nf -profile conda,slurm \
  --max_datasets 8 \
  --db_uri /scratch/multiomics/nickyoungblut/tiledb-loader/tiledb_prod3_tmp \
  --input_dir /processed_datasets/scRecount/scRecounter/prod3
```

```bash
nextflow run main.nf -profile conda,slurm \
  --max_datasets 8 \
  --db_uri /scratch/multiomics/nickyoungblut/tiledb-loader/tiledb_prod3_tmp \
  --input_dir /processed_datasets/scRecount/cellxgene/counted_SRXs
```

## Test scale

```bash
rm -rf /scratch/multiomics/nickyoungblut/tiledb-loader/tiledb_prod3_tmp
```

```bash
nextflow run main.nf -profile conda,slurm,report,trace \
  --max_datasets 100 \
  --mtx_batch_size 4 \
  --h5ad_batch_size 4 \
  --db_uri /scratch/multiomics/nickyoungblut/tiledb-loader/tiledb_tmp_4-4 \
  --input_dir /processed_datasets/scRecount/scRecounter/prod3
```

Time: 10m 15s

```bash
nextflow run main.nf -profile conda,slurm,report,trace \
  --max_datasets 100 \
  --mtx_batch_size 20 \
  --h5ad_batch_size 4 \
  --db_uri /scratch/multiomics/nickyoungblut/tiledb-loader/tiledb_tmp_20-4 \
  --input_dir /processed_datasets/scRecount/scRecounter/prod3
```

Time: 9m 55s (from-memory)
Time: 11m (from-disk)

```bash
nextflow run main.nf -profile conda,slurm,report,trace \
  --max_datasets 100 \
  --mtx_batch_size 50 \
  --h5ad_batch_size 4 \
  --db_uri /scratch/multiomics/nickyoungblut/tiledb-loader/tiledb_tmp_50-4 \
  --input_dir /processed_datasets/scRecount/scRecounter/prod3
```

Time: 10m 25s


```bash
nextflow run main.nf -profile conda,slurm,report,trace \
  --max_datasets 100 \
  --mtx_batch_size 25 \
  --h5ad_batch_size 2 \
  --db_uri /scratch/multiomics/nickyoungblut/tiledb-loader/tiledb_tmp_25-2 \
  --input_dir /processed_datasets/scRecount/scRecounter/prod3
```

Time: `TODO`