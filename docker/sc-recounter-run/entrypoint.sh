#!/bin/bash

# create run name
RUN_NAME="SCRECOUNTER_$(date +"%Y-%m-%d_%H-%M-%S")"

# Activate the micromamba environment and run the pipeline
micromamba run -n sc-recounter-run \
  nextflow run main.nf \
    -profile docker,gcp,report,trace \
    -name $RUN_NAME \
    -work-dir "gs://arc-ctc-nextflow/scRecounter/prod/work/${RUN_NAME}" \
    --output_dir "gs://arc-ctc-screcounter/prod/${RUN_NAME}" \
    -ansi-log false "$@"

# Remove the work directory
micromamba run -n sc-recounter-run \
  nextflow clean $RUN_NAME -force -quiet 
