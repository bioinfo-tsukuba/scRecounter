#!/bin/bash

# create run name
RUN_NAME="SCRECOUNTER_$(date +"%Y-%m-%d_%H-%M-%S")"

# Activate the micromamba environment and run the pipeline
micromamba run -n sc-recounter-run \
  nextflow run main.nf \
    -profile docker,gcp,prod,report,trace \
    -name $RUN_NAME \
    -work-dir gs://arc-ctc-nextflow/scRecounter/prod/work/${RUN_NAME} \
    -ansi-log false "$@"