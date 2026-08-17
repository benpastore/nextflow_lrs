FROM ubuntu:20.04

LABEL maintainer="https://github.com/twolinin/longphase"
LABEL version="2.0.2"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        git \
        g++ \
        gcc \
        autoconf \
        automake \
        make \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        ca-certificates \
        samtools \
        tabix \
        bcftools \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/longphase

RUN git clone https://github.com/twolinin/longphase.git /opt/longphase && \
    autoreconf -i && \
    ./configure && \
    make -j 4 && \
    rm -rf /opt/longphase/.git

ENV PATH="/opt/longphase:${PATH}"

CMD ["longphase", "phase", "--help"]