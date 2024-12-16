# Use micromamba base image
FROM mambaorg/micromamba:1.5.7

# Use bash shell
SHELL ["/bin/bash", "-c"]

# Set user to root for installation (already root by default, but kept for clarity)
USER root

# Set base directory of the Dockerfile; assuming build is run from the root of the repository
ARG BASE_DIR=docker/auto-run

# Install OS-level packages
RUN apt-get update -y \
    && apt-get install -y build-essential procps curl git \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*

# Install cloud_sql_proxy
# ARG PROXY_URL="https://storage.googleapis.com/cloud-sql-connectors/cloud-sql-proxy/v2.14.2/cloud-sql-proxy.linux.amd64"
# RUN curl $PROXY_URL -o /usr/bin/cloud_sql_proxy \
#     && chmod +x /usr/bin/cloud_sql_proxy \
#     && mkdir -p /cloudsql \
#     && chmod 777 /cloudsql

# Copy environment file
COPY --chown=$MAMBA_USER:$MAMBA_USER ${BASE_DIR}/environment.yml /tmp/environment.yml

# Install the environment using micromamba
ARG GITHUB_PAT
RUN micromamba create -n auto-run -f /tmp/environment.yml \
    && micromamba clean --all --yes \
    && micromamba run -n auto-run pip install git+https://${GITHUB_PAT}@github.com/ArcInstitute/SRAgent.git 

# Install the auto-run script
COPY --chown=$MAMBA_USER:$MAMBA_USER ${BASE_DIR}/auto-run.py /usr/bin/auto-run.py

# Copy the nextflow pipeline
#RUN mkdir -p /app
#COPY --chown=$MAMBA_USER:$MAMBA_USER main.nf nextflow.config workflows/ lib/ config/ bin/ /app/

# Copy the entrypoint script
#COPY ${BASE_DIR}/entrypoint.sh /usr/bin/entrypoint.sh
#RUN chmod +x /usr/bin/entrypoint.sh

# Activate the environment by default
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV ENV_NAME=auto-run

# Set the entrypoint script
CMD ["auto-run.py"]
