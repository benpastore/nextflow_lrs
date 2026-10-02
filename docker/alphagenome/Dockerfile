# Local (non-API) AlphaGenome inference, built from google-deepmind's
# research/local-inference repo -- NOT the google-deepmind/alphagenome pip
# package, which is a client for DeepMind's hosted API only (no local
# weights support, needs an API key + network access at inference time,
# which most HPC compute nodes don't have). This is what
# bin/run_alphagenome.py actually calls (dna_model.create_from_kaggle(...)
# against a pre-downloaded local weights cache -- see params.alphagenome_weights
# in workflows/haplotype_aware/nextflow.config).
#
# Install steps follow alphagenome_research's own README exactly
# (https://github.com/google-deepmind/alphagenome_research):
#   git clone https://github.com/google-deepmind/alphagenome_research.git
#   pip install -e ./alphagenome_research
# ...with one addition: that install pulls in a CPU-only `jax` wheel by
# default, so GPU support is installed explicitly after, per JAX's own
# install docs (https://docs.jax.dev/en/latest/installation.html):
#   pip install -U "jax[cuda12]"
# (the "jax[cuda12]" pip extra bundles its own CUDA/cuDNN via wheels, so a
# `-runtime` base image is sufficient here -- nothing needs to compile
# against the CUDA toolkit headers a `-devel` image would add.)
#
# requires-python >=3.11 per alphagenome_research's pyproject.toml --
# Ubuntu 24.04 ships Python 3.12 by default, so no deadsnakes/build-from-
# source needed.
#
# Build on OSC (no Docker daemon there -- buildah is OSC's documented
# daemonless, rootless-capable Dockerfile builder), then convert straight to
# a local .sif -- no registry/Docker Hub push needed:
#   buildah build -f docker/alphagenome/Dockerfile --format docker -t alphagenome:latest docker/alphagenome/
#   buildah push alphagenome:latest docker-archive:alphagenome.tar
#   apptainer build alphagenome.sif docker-archive://alphagenome.tar
# Then place alphagenome.sif at params.alphagenome_sif (see
# workflows/haplotype_aware/nextflow.config), same convention as
# params.herro_inference_sif.

FROM nvidia/cuda:12.6.3-cudnn-runtime-ubuntu24.04

LABEL maintainer="https://github.com/google-deepmind/alphagenome_research"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 https://github.com/google-deepmind/alphagenome_research.git /opt/alphagenome_research

WORKDIR /opt/alphagenome_research

# --break-system-packages: Ubuntu 24.04's system pip enforces PEP 668
# (externally-managed environment) by default -- fine to override in a
# single-purpose container with no other Python-using system packages.
RUN pip install --no-cache-dir --break-system-packages -e . \
    && pip install --no-cache-dir --break-system-packages -U "jax[cuda12]"

WORKDIR /

CMD ["python3", "-c", "from alphagenome_research.model import dna_model; print('alphagenome_research OK')"]
