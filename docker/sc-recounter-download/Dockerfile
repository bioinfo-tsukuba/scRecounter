# Use micromamba base image
FROM mambaorg/micromamba:1.5.7

# Use bash shell
SHELL ["/bin/bash", "-c"]

# Set user to root for installation (already root by default, but kept for clarity)
USER root

# Install OS-level packages
RUN apt-get update -y \
    && apt-get install -y build-essential procps curl \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*

# Copy environment file
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Install the environment using micromamba
RUN micromamba create -f /tmp/environment.yml \
    && micromamba clean --all --yes

# Activate the environment by default
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV ENV_NAME=download