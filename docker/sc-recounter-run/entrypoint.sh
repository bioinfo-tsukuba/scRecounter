#!/bin/bash

# Activate the micromamba environment and run the pipeline
micromamba run -n sc-recounter-run \
  nextflow run main.nf \
    -profile docker,gcp,dev,no_acc_dev,report,trace -ansi-log false "$@"