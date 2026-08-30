FROM condaforge/mambaforge:latest

# paraphase 3.5.0 (the pre-Rust-rewrite Python release -- 4.0.0's Rust
# realignment step panics on every CHM13 read, see PARAPHASE label in
# nextflow.config) shells out to external samtools/minimap2 binaries but
# doesn't declare them as bioconda dependencies, so they're not present in
# the stock quay.io/biocontainers/paraphase:3.5.0 image. Bundle them here.
RUN mamba install -y -c bioconda -c conda-forge \
        paraphase=3.5.0 \
        samtools \
        minimap2 \
    && mamba clean -afy

WORKDIR /data
