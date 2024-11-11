scRecouter
==========

Nextflow pipeline to re-process public single-cell RNA-seq data.

# Workflow

* User provides either:
  * A table of samples & associated accessions
  * A table of samples & downloaded fastq files
  * Associated files required:
    * a table of input samples and associated fastq files
    * a table of barcodes to use for cell barcode and UMI identification
    * a table of STAR index files to use for mapping
  * The user can also provide a table of STAR parameters, if they have already been determined
* Pipeline:
  * Read file formatting and QC
    * includes assessing read length
  * Cell barcode, UMI, and strand identification
    * via STARsolo with a subset of reads
      * assess fraction of valid barcodes
  * The STAR parameter setting with the highest fraction of valid barcodes is used for the full STAR run
    * What to do if there is a tie?
  * STAR run with selected parameters
    * All final count tables are "published" to the output directory

# Installation

## Conda & mamba install

`mamba` is needed to run the pipeline. 
It is a faster version of `conda`. `mamba` can be installed via `conda`. 

To install both `conda` and `mamba`, 
see the [conda/mamba Notion docs](https://www.notion.so/arcinstitute/Conda-Mamba-8106bed9553d46cca1af4e10f486bec2).

## Nextflow install

It is easiest to install Nextflow using `mamba`:

```bash
mamba install -n nextflow_env -c bioconda nextflow
```

Make sure to activate the environment before running the pipeline:

```bash
mamba activate nextflow_env
```

## Pipeline install

### Add ssh key to GitHub

> This is only needed if you have not already added your ssh key to GitHub.

```bash
ssh-keygen -t ed25519 -C "your_email@example.com"
```

* change `your_email@example.com` to your Arc email

```bash
cat ~/.ssh/id_ed25519.pub
```

* copy the output
* GoTo: `GitHub => Settings > SSH and GPG keys > New SSH key`
* Paste the output into the key field
* Add a title (e.g., `Chimera`)
* Click `Add SSH key`

### Clone the repo

```bash
git clone git@github.com:ArcInstitute/scRecounter.git \
  && cd scRecounter
```

### Pipeline conda environments 

The pipeline uses conda environments to manage dependencies. 
Nextflow will automatically create the environments as long as `mamba` is installed.

**Note:** it can take a while to create the environments, even with `mamba`.


# Usage

## Input

Input can either be:

1. `Accessions table` => A table of accessions per sample. The read data will be downloaded from the SRA.
1. `Reads table` => A table of read files per sample. The read data has already been downloaded.

Note: in either case, reads will be merged by `sample` if there are multiple read files (accessions) per sample.

### Accessions table

Lists the samples and their associated SRA experiment accessions.

Example:

| Sample    | Accession    |
|-----------|--------------|
| sample1   | SRR13112659  |
| sample1   | SRR13112660  |
| sample2   | SRR13112661  |


### Reads table

Lists the samples and their associated read (fastq) files.

Example:


| sample                  | fastq_1                                          | fastq_2                                          |
|-------------------------|--------------------------------------------------|--------------------------------------------------|
| SRX10188997_SRR13806043 | path/to/reads/SRX10188997_SRR13806043_1.fastq.gz | path/to/reads/SRX10188997_SRR13806043_2.fastq.gz |
| SRX10188963_SRR13806077 | path/to/reads/SRX10188963_SRR13806077_1.fastq.gz | path/to/reads/SRX10188963_SRR13806077_2.fastq.gz |


### Barcode table

Lists all of the possible barcodes that will be used to determine the cell barcode and UMI for the samples.

Example:

| name             | cell_barcode_length | umi_length | file_path                                                                |
|------------------|---------------------|------------|--------------------------------------------------------------------------|
| 737K-arc-v1      | 16                  | 12         | /large_storage/goodarzilab/public/scRecount/genomes/737K-arc-v1.txt      |
| 737K-august-2016 | 16                  | 12         | /large_storage/goodarzilab/public/scRecount/genomes/737K-august-2016.txt |
| 3M-february-2018 | 16                  | 10         | /large_storage/goodarzilab/public/scRecount/genomes/3M-february-2018.txt |


### STAR index table

Lists the STAR index files that will be used to map the reads.

Example:

| Organism      | Star Index Path                                                                   |
|---------------|-----------------------------------------------------------------------------------|
| Homo sapiens  | /large_storage/goodarzilab/public/scRecount/genomes/star_refData_2020_hg38        |
| Mus musculus  | /large_storage/goodarzilab/public/scRecount/genomes/star2.7.11_refData_2020_mm10  |


### STAR params table 

Useful if you have already determined the STAR parameters for each sample. 
Provide via `--star_params`. 

Example:

| Sample    | Fastq_1              | Fastq_2              | Barcodes File         | Star Index                   | Cell Barcode Length | UMI Length | Strand |
|-----------|----------------------|----------------------|-----------------------|------------------------------|---------------------|------------|--------|
| Sample_1  | sample_1_R1.fastq.gz | sample_1_R2.fastq.gz | sample_1_barcodes.tsv | star_refData_2020_hg38       | 16                  | 10         | +      |
| Sample_2  | sample_2_R1.fastq.gz | sample_2_R2.fastq.gz | sample_2_barcodes.tsv | star2.7.11_refData_2020_mm10 | 16                  | 10         | -      |

## Nextflow run 

### Characterize datasets

Use just a small subset of reads in the dataset to identify library prep method, species, etc.

TODO


```bash
nextflow run main.nf -profile conda,vm \
  --keep_temp true \
  --define \
  --max_spots 50000 \
  --accessions data/accessions_var.csv \
  --outdir tmp/accessions_var \
  -resume
```


# Dev

## Run

#### Local, defining STAR params

```bash
nextflow run main.nf \
  -profile conda,vm  \
  --keep_temp true  \
  --define \
  --max_spots 50000 \
  --accessions data/accessions_var.csv \
  --outdir tmp/accessions_var \
  -resume
```

Problematic accessions

```bash
nextflow run main.nf \
  -profile conda,vm  \
  --keep_temp true  \
  --define \
  --max_spots 50000 \
  --accessions data/accessions_problems.csv \
  --outdir tmp/accessions_problems \
  -resume
```

#### Local, with accessions

```bash
nextflow run main.nf -profile conda,vm,dev_acc
```

#### SLURM, with accessions

A couple of small-data accessions

```bash
nextflow run main.nf \
  -profile conda,slurm \
  --keep_temp true \
  --accessions data/accessions_small_n2.csv \
  --outdir tmp/results_small_n2
```

Many small-data accessions, subsampled

```bash
nextflow run main.nf \
  -profile conda,slurm \
  --keep_temp true \
  --max_spots 500000 \
  --accessions data/accessions_small_n10.csv \
  --outdir tmp/results_small_n10
```

## Convert accessions

Use the `acc2srr.py` script in this repo. An example:

```bash
./scripts/acc2srr.py --email YOUR_EMAIL@arcinstitute.org accessions.txt 
```

The `accessions.txt` file should contain a list of GEO or SRA accessions (one per line). 

## Convert Docker container to Apptainer

Pull the docker image (e.g., `ubuntu:latest`) and convert it to an Apptainer container:

```bash
apptainer pull ubuntu_latest.sif docker://ubuntu:latest
```

## GCP

Create VM

```bash
gcloud compute instances create sc-recounter-vm \
    --project=c-tc-429521 \
    --zone=us-east1-b \
    --machine-type=e2-standard-4 \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud \
    --boot-disk-size=50GB \
    --boot-disk-type=pd-balanced \
    --tags=allow-http,allow-https \
    --scopes=storage-full
```

ssh onto the VM

```bash
gcloud compute ssh sc-recounter-vm \
  --zone=us-east1-b \
  --project=c-tc-429521 \
  --impersonate-service-account=${HOME}/.gcp/c-tc-429521-6f6f5b8ccd93.json
```

Install micromamba

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

Install nextflow

```bash
mamba install -n base -c conda-forge -c bioconda nextflow
```

Run the pipeline

```bash
nextflow run main.nf -profile docker,gcp,dev_acc_gcp
```

To stop the VM 

```bash
gcloud compute instances stop sc-recounter-vm --zone=us-east1-b
```

## Resources

### Software

* [sra-tools](https://github.com/ncbi/sra-tools)
  * download data from the SRA
* [gget](https://github.com/pachterlab/gget)
  * efficient querying of genomic databases
* [ffq](https://github.com/pachterlab/ffq)
  * Fetch metadata information from the SRA and other databases
  * Can be used to get SRA study accessions from paper DOIs
* [pysradb](https://github.com/saketkc/pysradb)
  * A Python package for retrieving metadata from SRA/ENA/GEO
* [gencube](https://github.com/snu-cdrc/gencube)
  * Efficient retrieval, download, and unification of genomic data from leading biodiversity databases
* [geofetch](https://pep.databio.org/geofetch/)
  * Downloads and organizes data and metadata from GEO and SRA
* [nf-core/fetchngs](https://nf-co.re/fetchngs/1.12.0/)
  * Nextflow pipeline for downloading NGS data

### Databases

* [Single cell studies database](https://docs.google.com/spreadsheets/d/1En7-UV0k0laDiIfjFkdn7dggyR7jIk3WH8QgXaMOZF0/edit?gid=0#gid=0)

## Workflow

* Input
  * csv of SRA experiment accessions
* Download
  * adapted from [nf-core/fetchngs](https://nf-co.re/fetchngs/1.12.0/)
    * https://github.com/nf-core/fetchngs/blob/master/workflows/sra/main.nf
* Read QC 
  * seqkit stats
* Characterize datasets
  * See https://github.com/ArcInstitute/scRecount/blob/chris/chris_scripts/dsub_solution/process_srr.sh
* Map to reference
  * STARsolo
  * See https://github.com/ArcInstitute/scRecount/blob/chris/chris_scripts/dsub_solution/process_srr.sh