scRecounter
===========

# This code folked from https://github.com/ArcInstitute/scRecounter

A Nextflow pipeline to re-process single-cell RNA-seq data from the Sequence Read Archive (SRA) with integrated process tracking and error handling.

# Workflow

* **User provides:**
  * A table of samples & associated accessions (CSV format)
  * Associated files required:
    * A table of barcodes to use for cell barcode and UMI identification (`data/barcodes.csv`)
    * A table of STAR index directories to use for mapping (`data/star_indices.csv`)
    * A .env.local file including DB information

* **Pipeline:**
  * Initialize process tracking for experiment monitoring in SQL
  * Load accessions from provided table
  * Run `sra-stat` to check file sizes and filter out files exceeding the maximum size limit
  * For each accession:
    * Use `fastq-dump` to download a subset of reads as fastq files from the SRA
    * Determine the "best" STAR parameters by mapping the reads using various parameter combinations
      * Parameters: version of cell barcodes, cell barcode length, UMI length, strand, STAR reference index
      * The STAR parameters are selected based on the fraction of valid barcodes
    * Download all reads with `fasterq-dump`
      * If download fails, try again with `fastq-dump` using a max of `fallback_max_spots` reads (see `nextflow.config`)
    * Map the reads with STARsolo using the "best" STAR parameters
  * Track individual accession results and completion status in SQL

# Manuscript
## Original paper
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
git clone https://github.com/bioinfo-tsukuba/scRecounter \
  && cd scRecounter
```

### Pipeline conda environments (if running locally)

The pipeline uses conda environments to manage dependencies. 
Nextflow will automatically create the environments as long as `mamba` is installed.

**Note:** it can take a while to create the environments, even with `mamba`.

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

### Accessions table (Required)

CSV file listing samples and their associated SRA experiment accessions.

> Not required if pulling accessions from the scRecounter SQL database.
  Leave `--accessions` empty or omit it to use database accessions.

Example format:

| sample      | accession   | organism |
|-------------|-------------|----------|
| SRX22716300 | SRR27024456 | human    |
| SRX25994842 | SRR30571763 | mouse    |

or 

| sample | accession | download_url | organism |
|-------------|-------------|----------|----------|
| SRX12280794 | SRR15992285 | https://ddbj.nig.ac.jp/public/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX122/SRX12280794/SRR15992285/SRR15992285.sra | human
| SRX12101437 | SRR15808974 | https://ddbj.nig.ac.jp/public/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX121/SRX12101437/SRR15808974/SRR15808974.sra | human
| SRX12101437 | SRR15808975 | https://ddbj.nig.ac.jp/public/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX121/SRX12101437/SRR15808975/SRR15808975.sra | human

If you provide download_url, scRecounter tries to download from the url.

### Barcode table (Required)

CSV file (`data/barcodes.csv`) listing cell barcode files for different 10x chemistry versions.

Example format:

| name             | cell_barcode_length | umi_length | file_path                                                                |
|------------------|---------------------|------------|--------------------------------------------------------------------------|
| 737K-arc-v1      | 16                  | 12         | /path/to/barcodes/737K-arc-v1.txt      |
| 737K-august-2016 | 16                  | 12         | /path/to/barcodes/737K-august-2016.txt |
| 3M-february-2018 | 16                  | 10         | /path/to/barcodes/3M-february-2018.txt |

### STAR index table (Required)

CSV file (`data/star_indices.csv`) listing pre-built STAR genome indices.

Example format:

| Organism | Star Index Path                                                                   |
|----------|-----------------------------------------------------------------------------------|
| human    | /path/to/genomes/Index        |
| mouse    | /path/to/genomes/Index  |

You need to prepare index files using STAR. 
We downloaded the genome annotations from gencode.

Prepare STAR index
```bash
STAR --runMode genomeGenerate \
    --runThreadN 20 \
    --genomeDir "Index/" \
    --genomeFastaFiles "path/to/genome.fa" \
    --genomeSAindexNbases 14 \
    --genomeChrBinNbits 18 \
    --genomeSAsparseD 3 \
    --limitGenomeGenerateRAM 17179869184 \
    --sjdbGTFfile "path/to/gtf" \
    --sjdbOverhang 100
```

## Process Tracking

scRecounter includes PostgreSQL-based process tracking for monitoring pipeline execution status. This feature tracks experiment progress and error states in a local database.

### Process Tracking Setup (Required)

1. **Prepare PostgreSQL:**
```bash
docker pull postgres:16
docker compose up -d
```

2. **Create database and user:**
```bash
# Add account
psql "host=localhost user=admin dbname=appdb password=PLEASE_SET_ADMIN_PASSWARD" -c "CREATE USER YOURACCOUNT WITH PASSWORD 'YOURPASSWORD';"
# Add table
psql "host=localhost user=admin dbname=appdb password=PLEASE_SET_ADMIN_PASSWARD" -c "CREATE TABLE IF NOT EXISTS public.experiment_process (id text primary key, experiment_id text, srx_accession text, organism text, analysis_date date, process_type text, status integer, start_datetime timestamp, finish_datetime timestamp, path text, process_id text, accessions_file text);"
# Add permission
psql "host=localhost user=admin dbname=appdb password=PLEASE_SET_ADMIN_PASSWARD" -c "GRANT SELECT, INSERT, UPDATE, DELETE ON public.experiment_process TO YOURACCOUNT;"
```

3. **Configure environment:**
The `.env.local` file contains database connection settings:
```
LOCAL_DB_HOST=localhost
LOCAL_DB_NAME=appdb
LOCAL_DB_USER=USER_NAME
LOCAL_DB_PASSWORD=PASSWORD
LOCAL_DB_PORT=PORT_ID
```

4. **Enable tracking:**
Process tracking is currently implemented but commented out in `main.nf`. To enable, uncomment the relevant PROCESS_TRACKER sections.

For detailed process tracking documentation, see [PROCESS_TRACKER_INTEGRATION.md](./PROCESS_TRACKER_INTEGRATION.md).

## Running the Pipeline

### Basic Usage

**Local execution with provided accessions:**
### Run command example in Cell-IO mapping
**Small test:**
Use sample ID
```bash
nextflow run main.nf \
  -work-dir tmp/work \
  -profile conda,trace,report,vm,vm_dev,dev,acc_dev_problems
  --barcodes     path/to/barcodes.csv     \
  --star_indices path/to/star_indices.csv \
  --accessions   data/accessions_small_n3.csv \
  --output_dir   OUTPUTDIR
```

Use URL
```bash
nextflow run main.nf \
  -work-dir tmp/work \
  -profile conda,trace,report,vm,vm_dev,dev,acc_dev_problems
  --barcodes     path/to/barcodes.csv     \
  --star_indices path/to/star_indices.csv \
  --accessions   data/accessions_url_n6.csv \
  --output_dir   OUTPUTDIR
```

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
