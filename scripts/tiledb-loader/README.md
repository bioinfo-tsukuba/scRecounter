tiledb loader
=============

# Production

prod3

```bash
nextflow run main.nf -profile conda,slurm,dev \
  --db_uri  /scratch/multiomics/nickyoungblut/tiledb-loader/tiledb-soma_prod3_GeneFull_Ex50pAS \
  --input_dir /processed_datasets/scRecount/scRecounter/prod3
```

```bash
nextflow run main.nf -profile conda,slurm,dev \
  --db_uri  /scratch/multiomics/nickyoungblut/tiledb-loader/tiledb-soma_prod3_GeneFull_Ex50pAS \
  --input_dir /processed_datasets/scRecount/cellxgene/counted_SRXs
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
nextflow run main.nf -profile conda,slurm \
  --max_datasets 100 \
  --mtx_batch_size 8 \
  --h5ad_batch_size 2 \
  --db_uri /scratch/multiomics/nickyoungblut/tiledb-loader/tiledb_prod3_tmp \
  --input_dir /processed_datasets/scRecount/scRecounter/prod3
```

Time: 

```bash
nextflow run main.nf -profile conda,slurm \
  --max_datasets 100 \
  --mtx_batch_size 8 \
  --h5ad_batch_size 4 \
  --db_uri /scratch/multiomics/nickyoungblut/tiledb-loader/tiledb_prod3_tmp \
  --input_dir /processed_datasets/scRecount/scRecounter/prod3
```

Time: 