FROM staphb/minimap2:latest

USER root

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        samtools \
        bcftools \
        tabix \
        ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

CMD ["minimap2"]