FROM condaforge/mambaforge:latest

# lbcb-sci/herro (https://github.com/lbcb-sci/herro) -- CPU-side steps only
# (preprocess.sh / create_batched_alignments.sh). GPU inference itself uses
# the maintainers' own prebuilt herro.sif (see nextflow.config, label
# 'herro_inference') rather than being built here.
RUN apt-get update && apt-get install -y --no-install-recommends git build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 https://github.com/lbcb-sci/herro.git /opt/herro

# exact env spec from the repo (scripts/herro-env.yml) -- includes
# seqkit/minimap2 (conda) and duplex-tools + a patched porechop fork (pip),
# both of which preprocess.sh shells out to.
RUN mamba env create -n herro -f /opt/herro/scripts/herro-env.yml \
    && mamba clean -afy

ENV PATH=/opt/conda/envs/herro/bin:$PATH

WORKDIR /data
