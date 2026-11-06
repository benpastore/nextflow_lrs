process STRAGLR {

    tag "$sampleID"
    label 'straglr'

    publishDir "${params.results}/06_variants/straglr", mode: 'copy'

    input:
        tuple val(sampleID), val(bam), val(bai)
        tuple val(ref_fa), val(ref_fai)
        
    output:
        tuple val(sampleID),
            path("*.straglr.tsv"),
            path("*.straglr.bed"),
            emit: straglr_ch

    script:
    def loci_arg = params.straglr_loci_bed ? "--loci ${params.straglr_loci_bed}" : ""

    """
    #!/bin/bash
    set -euo pipefail

    name=\$(basename ${bam} .bam)

    # straglr's tre.py encodes internal state as read_id:qstart:qend:tpos:tlen
    # and splits on ':', assuming read IDs never contain one -- but herro
    # names split-fragment reads like <uuid>:0/:1/... (see MINIMAP2_ALIGN's
    # -y comment for the same underlying quirk), which breaks that split
    # ("too many values to unpack"). No samtools binary in this container
    # (straglr's own conda deps are just trf/blast/pysam/pybedtools), so
    # sanitize QNAMEs with pysam directly rather than shelling out.
    python3 - <<PYEOF
import pysam
inp = pysam.AlignmentFile("${bam}", "rb")
out = pysam.AlignmentFile("\${name}.qname_sanitized.bam", "wb", template=inp)
for read in inp:
    if ":" in read.query_name:
        read.query_name = read.query_name.replace(":", "_")
    out.write(read)
inp.close()
out.close()
pysam.index("\${name}.qname_sanitized.bam")
PYEOF

    straglr.py \\
        \${name}.qname_sanitized.bam \\
        ${ref_fa} \\
        \${name}.straglr \\
        ${loci_arg} \\
        --genotype_in_size \\
        --min_support ${params.straglr_min_support ?: 2} \\
        --max_str_len ${params.straglr_max_str_len ?: 100} \\
        --max_num_clusters ${params.straglr_max_num_clusters ?: 2} \\
        --nprocs ${task.cpus} \\
        ${params.straglr_args ?: ''}

    # sanity check: straglr's BED must start with its own "#chrom" header --
    # a run once produced a BED whose header had a filename glued in front
    # of it ("<prefix>.straglr.haplotagged.tsv#chrom..."); fail here, in the
    # task that wrote it, rather than letting it reach downstream parsers
    for f in \${name}.straglr.bed \${name}.straglr.tsv; do
        ls -l --time-style=full-iso "\$f" >&2
    done
    first=\$(head -n 1 \${name}.straglr.bed)
    if [[ "\$first" != "#chrom"* ]]; then
        echo "[STRAGLR] \${name}.straglr.bed header is malformed: \${first:0:200}" >&2
        exit 1
    fi
    """
}