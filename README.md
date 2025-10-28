scRecounter
===========

A Nextflow pipeline to re-process single-cell RNA-seq data from the Sequence Read Archive (SRA) with integrated process tracking and error handling.

# Workflow

* **User provides:**
  * A table of samples & associated accessions (CSV format)
    * Alternatively, the pipeline can pull accessions from the scRecounter SQL database
  * Associated files required:
    * A table of barcodes to use for cell barcode and UMI identification (`barcodes.csv`)
    * A table of STAR index directories to use for mapping (`star_indices.csv`)
* **Pipeline:**
  * Initialize process tracking for experiment monitoring (optional)
  * Load accessions from provided table or SQL database
  * Run `sra-stat` to check file sizes and filter out files exceeding the maximum size limit
  * For each accession:
    * Use `fastq-dump` to download a subset of reads as fastq files from the SRA
    * Determine the "best" STAR parameters by mapping the reads using various parameter combinations
      * Parameters: version of cell barcodes, cell barcode length, UMI length, strand, STAR reference index
      * The STAR parameters are selected based on the fraction of valid barcodes
    * Download all reads with `fasterq-dump`
      * If download fails, try again with `fastq-dump` using a max of `fallback_max_spots` reads (see `nextflow.config`)
    * Map the reads with STARsolo using the "best" STAR parameters
  * Track individual accession results and completion status

# Manuscript

**scBaseCamp: An AI agent-curated, uniformly processed, and continually expanding single cell data repository**.
Nicholas D Youngblut, Christopher Carpenter, Jaanak Prashar, Chiara Ricci-Tam, Rajesh Ilango, Noam Teyssier,
Silvana Konermann, Patrick Hsu, Alexander Dobin, David P Burke, Hani Goodarzi, Yusuf H Roohani.
bioRxiv 2025.02.27.640494; doi: [https://doi.org/10.1101/2025.02.27.640494](https://doi.org/10.1101/2025.02.27.640494)

# Installation

## Conda & mamba install

`mamba` is needed to run the pipeline. 
It is a faster version of `conda`. 
`mamba` can be installed via `conda`. 
You can use `conda` instead of `mamba` if you prefer.

## Nextflow install

It is easiest to install Nextflow using `mamba` (or `conda`).

```bash
mamba create -n nextflow_env -c bioconda nextflow
```

Make sure to activate the environment before running the pipeline:

```bash
mamba activate nextflow_env
```

All other dependencies will be installed by Nextflow.


## Pipeline install

### Clone the repo

```bash
git clone https://github.com/ArcInstitute/scRecounter.git \
  && cd scRecounter
```

### Pipeline conda environments (if running locally)

The pipeline uses conda environments to manage dependencies. 
Nextflow will automatically create the environments as long as `mamba` is installed.

**Note:** it can take a while to create the environments, even with `mamba`.

### Pipeline Docker containers (if running on GCP) 

The pipeline defaults to using custom Docker containers hosted on Google Artifact Registry.

You can build the Docker containers yourself. See [./docker/README.md](./docker/README.md) for details.
Be sure to update the [profiles.config](./config/profiles.config) file to point to the new containers.

# Usage

## Configuration Parameters

Key parameters in `nextflow.config`:
- `max_samples`: Maximum number of samples to process (default: 3)
- `max_accessions`: Maximum accessions per sample for STAR parameter testing (default: 1) 
- `max_spots`: Maximum reads for STAR parameter assessment (default: 1,000,000)
- `max_sra_size`: Maximum SRA file size in GB to process (default: 300GB)
- `min_read_len`: Minimum read length threshold (default: 26bp)
- `organisms`: Target organisms when pulling from database (default: "human,mouse")

## Input Files

### Accessions table (Optional)

CSV file listing samples and their associated SRA experiment accessions.

> Not required if pulling accessions from the scRecounter SQL database.
  Leave `--accessions` empty or omit it to use database accessions.

Example format:

| sample      | accession   | organism |
|-------------|-------------|----------|
| SRX22716300 | SRR27024456 | human    |
| SRX25994842 | SRR30571763 | mouse    |

> The `organism` column is optional and helps reduce parameter search space by pre-filtering STAR indices.

### Barcode table (Required)

CSV file (`data/barcodes.csv`) listing cell barcode files for different 10x chemistry versions.

Example format:

| name             | cell_barcode_length | umi_length | file_path                                                                |
|------------------|---------------------|------------|--------------------------------------------------------------------------|
| 737K-arc-v1      | 16                  | 12         | /large_storage/goodarzilab/public/scRecount/genomes/737K-arc-v1.txt      |
| 737K-august-2016 | 16                  | 12         | /large_storage/goodarzilab/public/scRecount/genomes/737K-august-2016.txt |
| 3M-february-2018 | 16                  | 10         | /large_storage/goodarzilab/public/scRecount/genomes/3M-february-2018.txt |

### STAR index table (Required)

CSV file (`data/star_indices.csv`) listing pre-built STAR genome indices.

Example format:

| Organism | Star Index Path                                                                   |
|----------|-----------------------------------------------------------------------------------|
| human    | /large_storage/goodarzilab/public/scRecount/genomes/star_refData_2020_hg38        |
| mouse    | /large_storage/goodarzilab/public/scRecount/genomes/star2.7.11_refData_2020_mm10  |

> When `organism` is specified in the accessions table, only matching STAR indices are tested, reducing computation time.

## Running the Pipeline

### Basic Usage

**Local execution with provided accessions:**
```bash
nextflow run main.nf \
  -work-dir tmp/work \
  -profile conda,trace,report,vm,vm_dev,dev,acc_dev \
  --accessions data/my_accessions.csv
```

**Local execution using database accessions:**
```bash
nextflow run main.nf \
  -work-dir tmp/work \
  -profile conda,trace,report,vm,vm_dev,dev,no_acc_dev
```

**GCP execution with Docker containers:**
```bash
nextflow run main.nf \
  -profile docker,trace,report,gcp,gcp_dev,dev,acc_dev \
  --accessions gs://my-bucket/accessions.csv \
  --output_dir gs://my-bucket/results/
```

### Testing and Development

**Small test with problematic datasets:**
```bash
nextflow run main.nf \
  -work-dir tmp/work \
  -profile conda,trace,report,vm,vm_dev,dev,acc_dev_problems
```

**Parameter characterization (small subset):**
```bash
nextflow run main.nf \
  -work-dir tmp/work \
  -profile conda,trace,report,vm,vm_dev,dev \
  --max_spots 100000 \
  --accessions data/test_accessions.csv \
  --define true
```

### Advanced Configuration

**Custom parameters:**
```bash
nextflow run main.nf \
  -profile conda,trace,report \
  --max_samples 10 \
  --max_sra_size 500 \
  --min_read_len 30 \
  --organisms "human" \
  --output_dir custom_results/
```

### Cloud Deployment

For GCP Cloud Run deployment, see [./docker/sc-recounter-run/README.md](./docker/sc-recounter-run/README.md).


# Process Tracking (Optional)

scRecounter includes optional PostgreSQL-based process tracking for monitoring pipeline execution status. This feature tracks experiment progress and error states in a local database.

## Process Tracking Setup

1. **Install PostgreSQL:**
```bash
sudo apt update
sudo apt install postgresql postgresql-contrib
```

2. **Create database and user:**
```bash
sudo -u postgres createdb experimentprocess
sudo -u postgres createuser cellio
sudo -u postgres psql -c "ALTER USER cellio PASSWORD 'cEllIo_process';"
sudo -u postgres psql -c "GRANT ALL PRIVILEGES ON DATABASE experimentprocess TO cellio;"
```

3. **Configure environment:**
The `.env.local` file contains database connection settings:
```
LOCAL_DB_HOST=localhost
LOCAL_DB_NAME=experimentprocess
LOCAL_DB_USER=cellio
LOCAL_DB_PASSWORD=cEllIo_process
LOCAL_DB_PORT=5432
```

4. **Enable tracking:**
Process tracking is currently implemented but commented out in `main.nf`. To enable, uncomment the relevant PROCESS_TRACKER sections.

For detailed process tracking documentation, see [PROCESS_TRACKER_INTEGRATION.md](./PROCESS_TRACKER_INTEGRATION.md).

# Output Structure

Results are organized in the `output_dir` (default: `results/`):

```
results/
├── STAR/
│   └── [SAMPLE_ID]/
│       ├── Gene/
│       ├── GeneFull/
│       ├── GeneFull_Ex50pAS/
│       ├── GeneFull_ExonOverIntron/
│       ├── Velocyto/
│       └── [ACCESSION_ID]/
│           ├── merged_star_params.csv
│           └── selected_star_params.json
├── logs/
│   └── [SAMPLE_ID]/
│       └── [workflow_logs].log
├── nf-report/
└── nf-trace/
```

# Contributing

Feel free to fork the repository and submit a pull request. 
The priority is maintaining compatibility with the ongoing scBaseCamp project.