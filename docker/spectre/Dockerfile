FROM python:3.10-slim

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/venv/bin:${PATH}"

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    wget \
    curl \
    ca-certificates \
    build-essential \
    samtools \
    tabix \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

############################################################
# Install mosdepth from official GitHub release
############################################################

ARG MOSDEPTH_VERSION=0.3.11

RUN wget -q \
    https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth \
    -O /usr/local/bin/mosdepth && \
    chmod 755 /usr/local/bin/mosdepth

############################################################
# Python environment
############################################################

RUN python -m venv /opt/venv && \
    pip install --upgrade pip setuptools wheel build

############################################################
# Spectre dependencies
############################################################

RUN pip install \
    numpy==1.26.4 \
    pandas \
    scipy \
    matplotlib \
    pysam

############################################################
# Install Spectre from GitHub
############################################################

WORKDIR /opt

RUN git clone https://github.com/fritzsedlazeck/Spectre.git

WORKDIR /opt/Spectre

RUN python -m build && \
    pip install dist/*.tar.gz

############################################################
# Validation
############################################################

RUN python - <<'PY'
import numpy as np
print("numpy:", np.__version__)
assert int(np.__version__.split('.')[0]) < 2
PY

RUN mosdepth --help >/dev/null
RUN spectre --help >/dev/null
RUN samtools --version >/dev/null
RUN tabix --version >/dev/null

ENTRYPOINT ["/bin/bash"]