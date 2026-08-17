FROM continuumio/miniconda3:latest

# HERO.py shells out to `less {file} | head -1` to sniff fastq vs fasta
# (see fq_or_fa() in its source) -- the miniconda3 base image doesn't ship
# `less`, which makes that popen() return empty and crashes HERO with an
# IndexError on the first file it touches.
RUN apt-get update && apt-get install -y --no-install-recommends less \
    && rm -rf /var/lib/apt/lists/*

# HERO (https://github.com/kangxiongbin/HERO) has no bioconda/pip package of its own --
# its dependencies do, so we conda-install those per the upstream README and just git
# clone the tool itself alongside them.
RUN conda create -y -n hero -c bioconda -c conda-forge \
        python=3.6 racon minimap2 ropebwt2 fmlrc2 ratatosk lordec \
    && conda clean -afy

RUN git clone --depth 1 https://github.com/kangxiongbin/HERO.git /opt/HERO

ENV PATH=/opt/conda/envs/hero/bin:$PATH
ENV HERO_BIN=/opt/HERO/bin

WORKDIR /data
