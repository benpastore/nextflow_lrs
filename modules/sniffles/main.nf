
process SNIFFLES {

    label 'sniffles'
    tag "${sampleID}"

    input:
        tuple val(sampleID), val(bam), val(bai)

    output:
        tuple val(sampleID), path("*.sniffles.vcf"), emit : sniffles_raw_ch

    script:
    """
    #!/bin/bash

    name=\$(basename ${bam} .bam)

    # --phase uses the HP/PS tags already on this BAM (haplotagged by
    # longphase upstream) to emit a per-SV PHASE= field in INFO.
    #
    # bgzip/tabix aren't in this container (the biocontainers sniffles
    # image only ships sniffles itself, no htslib CLI tools) -- compression
    # + indexing is done by INDEX_SNIFFLES_VCF (bcftools container) instead.
    #
    # Sniffles 2.3.2 hangs at shutdown (multiprocessing workers never
    # join) AFTER writing a complete VCF -- jobs sat until the 12h SLURM
    # limit and died with exit 140. So run it in the background, wait for
    # its "Wrote N called SVs" line, then kill it and validate the VCF.
    sniffles \
      --input ${bam} \
      --vcf \$name.sniffles.vcf \
      --threads ${task.cpus} \
      --phase \
      --allow-overwrite > sniffles.log 2>&1 &
    pid=\$!

    while kill -0 \$pid 2>/dev/null; do
        if grep -q 'Wrote [0-9]* called SVs' sniffles.log; then
            sleep 30
            kill \$pid 2>/dev/null || true
            break
        fi
        sleep 10
    done
    wait \$pid
    rc=\$?

    cat sniffles.log

    if grep -q 'Wrote [0-9]* called SVs' sniffles.log \
        && [[ -s \$name.sniffles.vcf ]] \
        && [[ "\$(tail -c 1 \$name.sniffles.vcf | od -An -c | tr -d ' ')" == '\\n' ]]; then
        exit 0
    fi
    echo "sniffles did not finish writing \$name.sniffles.vcf (exit \$rc)" >&2
    exit \$rc
    """
}

process INDEX_SNIFFLES_VCF {

    tag "${sampleID}"
    label "bcftools"
    publishDir "${params.results}/06_variants/sniffles", mode: params.publish_mode

    input:
        tuple val(sampleID), path(vcf)

    output:
        tuple val(sampleID),
            path("*.sniffles.vcf.gz"),
            path("*.sniffles.vcf.gz.tbi"), emit : sniffles_ch

    script:
    """
    #!/bin/bash

    set -euo pipefail

    if [[ "${vcf}" == *.vcf.gz ]]; then
        tabix -f -p vcf ${vcf}
    else
        bgzip -f ${vcf}
        tabix -f -p vcf ${vcf}.gz
    fi
    """
}