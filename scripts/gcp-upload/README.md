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

## prod

### SRA

```bash
nextflow run main.nf -profile conda,slurm
```

### CZI

```bash
nextflow run main.nf -profile conda,slurm --input_dir /processed_datasets/scRecount/cellxgene/counted_SRXs
```

