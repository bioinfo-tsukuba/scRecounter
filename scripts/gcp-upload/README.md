gcp-loader
==========

A simple Nextflow pipeline for efficiently loading single-cell data as h5ad files onto GCP


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
nextflow run main.nf -profile conda,vm \
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

