Development docs
================


## deploy new GCP Cloud Run batch revision

See [./docker/sc-recounter/README.md](./docker/sc-recounter/README.md) for details.

## set env variables

```bash
export GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp/c-tc-429521-6f6f5b8ccd93.json"
```

## Run locally

Basic run:

```bash
nextflow run main.nf \
  -profile conda,vm,dev,vm_dev,acc_dev \
  -resume 
```

Run with problematic accessions

```bash
nextflow run main.nf \
  -profile conda,vm  \
  --define \
  --max_spots 50000 \
  --accessions data/accessions_problems.csv \
  --output_dir tmp/accessions_problems \
  -resume
```

Run with all organisms

```bash
nextflow run main.nf \
  -profile conda,slurm,dev,acc_all_org \
  -resume
```


## SLURM, with defined accessions

A couple of small-data accessions

```bash
nextflow run main.nf \
  -profile conda,slurm \
  --accessions data/accessions_small_n2.csv \
  --output_dir tmp/results_small_n2
```

Many small-data accessions, subsampled

```bash
nextflow run main.nf \
  -profile conda,slurm \
  --max_spots 500000 \
  --accessions data/accessions_small_n10.csv \
  --output_dir tmp/results_small_n10
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

## GCP VM setup

Create VM

```bash
gcloud compute instances create sc-recounter-vm \
    --project=c-tc-429521 \
    --zone=us-east1-c \
    --machine-type=e2-standard-8 \
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
  --zone=us-east1-c \
  --project=c-tc-429521 \
  --impersonate-service-account=${HOME}/.gcp/c-tc-429521-6f6f5b8ccd93.json
```

Install micromamba

```bash
curl -L \
  -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" \
  && bash Miniforge3-$(uname)-$(uname -m).sh
```

Add the required channels

```bash
conda config --add channels nodefaults \
  && conda config --add channels pytorch \
  && conda config --add channels bioconda \
  && conda config --add channels conda-forge
```

Install nextflow

```bash
conda install -n nextflow-env nextflow \
  && conda activate nextflow-env
```

> Be sure that the `GOOGLE_APPLICATION_CREDENTIALS` env variable is set to the path of the service account key

The service account needs the following roles:

```bash
# batch
gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/batch.serviceAgent"

# compute
gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/compute.admin"

gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/compute.instanceAdmin.v1"

# storage
gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/storage.objectViewer"

# Artifact Registry
gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/artifactregistry.reader"

# network and Monitoring Access
gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/iam.serviceAccountUser"

gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/compute.networkUser"

gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/logging.viewer"

gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/monitoring.viewer"
```

Run the pipeline

```bash
nextflow run main.nf -profile docker,gcp,acc_dev
```

To stop the VM:

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

