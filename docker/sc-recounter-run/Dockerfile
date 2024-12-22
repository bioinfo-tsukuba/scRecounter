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
RUN micromamba create -n sc-recounter-run -f /tmp/environment.yml \
    && micromamba clean --all --yes \
    && rm -rf /opt/conda/pkgs/*

# Activate environment by default
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV ENV_NAME=sc-recounter-run

# Copy Nextflow pipeline and the Python runner script
COPY --chown=$MAMBA_USER:$MAMBA_USER main.nf nextflow.config workflows/ lib/ config/ bin/ .

# Copy entrypoint script
COPY bin/db_utils.py ${BASE_DIR}/sc-recounter-run.py ${BASE_DIR}/entrypoint.sh ./

# Create a directory for the mamba cache
RUN mkdir -p /.cache/mamba/ && chmod -R ugo+rwx /.cache/mamba/

# Set user to mamba
ENTRYPOINT ["bash", "entrypoint.sh"]