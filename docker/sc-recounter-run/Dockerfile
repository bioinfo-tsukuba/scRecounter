# Use micromamba base image
FROM mambaorg/micromamba:1.5.7

# Use bash shell
SHELL ["/bin/bash", "-c"]

# Set working directory
WORKDIR /app

# Set user to root for installation
USER root

# Install OS-level packages (if needed)
RUN apt-get update -y \
    && apt-get install -y build-essential procps curl \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*

# variables
ARG BASE_DIR="docker/sc-recounter-run"

# Copy environment file into container
COPY --chown=$MAMBA_USER:$MAMBA_USER ${BASE_DIR}/environment.yml /tmp/environment.yml

# Create conda environment with micromamba
RUN micromamba create -n sc-recounter-run -f /tmp/environment.yml --quiet \
    && micromamba clean --all --yes \
    && rm -rf /opt/conda/pkgs/*

# Activate environment by default
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV ENV_NAME=sc-recounter-run

# Copy Nextflow pipeline and the Python runner script
COPY --chown=$MAMBA_USER:$MAMBA_USER main.nf nextflow.config .
COPY --chown=$MAMBA_USER:$MAMBA_USER bin/ ./bin/
COPY --chown=$MAMBA_USER:$MAMBA_USER config/ ./config/
COPY --chown=$MAMBA_USER:$MAMBA_USER data/ ./data/
COPY --chown=$MAMBA_USER:$MAMBA_USER lib/ ./lib/ 
COPY --chown=$MAMBA_USER:$MAMBA_USER workflows/ ./workflows/

# Copy runner scripts
COPY bin/db_utils.py ${BASE_DIR}/entrypoint.sh ${BASE_DIR}/cleanup.py  ./

# Create a directory for the mamba cache
RUN mkdir -p /.cache/mamba/ /app/.nextflow/ /scratch/ \
  && chmod -R ugo+rwx /.cache/mamba/ /app/.nextflow/ /scratch/

# Set the NXF_HOME environment variable
ENV NXF_HOME=/app/.nextflow

# Set user to mamba
ENTRYPOINT ["bash", "entrypoint.sh"]