#!/bin/bash

# Activate the micromamba environment and forward arguments to the Python script
micromamba run -n sc-recounter-run python sc-recounter-run.py "$@"