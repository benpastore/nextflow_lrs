process CONSOLIDATE_VARIANTS {

    tag "$sampleID"
    label 'low'
    errorStrategy 'ignore'

    publishDir "${params.results}/06_variants/consolidated", mode: params.publish_mode

    // Trio mode: father_vcf / mother_vcf (parents' Clair3) and sv_origin
    // (bin/sv_trio_origin.py) add INFO/ORIGIN -- the parental haplotype each
    // variant is on -- to every site, plus <sample>.phase_block_origin.tsv.
    // Without parents they are NO_FILE and the master VCF has no ORIGIN.
    input:
        tuple val(sampleID),
            val(clair3_longphase_vcf),
            val(sniffles_vcf),
            val(dipcall_vcf),
            val(hapdiff_vcf),
            val(father_vcf),
            val(mother_vcf),
            val(sv_origin)
        path(script)

    output:
        tuple val(sampleID), path("*.master.vcf"), emit: consolidated_ch
        tuple val(sampleID), path("*.phase_block_origin.tsv"), optional: true, emit: phase_block_origin_ch

    script:
    def trio = file(father_vcf.toString()).name != 'NO_FILE' && file(mother_vcf.toString()).name != 'NO_FILE'
    def trio_args = trio ? "--father-vcf ${father_vcf} --mother-vcf ${mother_vcf} --sv-origin ${sv_origin} --block-output ${sampleID}.phase_block_origin.tsv" : ""
    """
    #!/bin/bash
    set -euo pipefail

    # conda's own activate.d hooks for this env (specifically
    # binutils_linux-64's, which references \$ADDR2LINE with no default)
    # aren't written to survive `set -u` -- drop it just for activation.
    set +u
    source activate rnaseq
    set -u

    python3 ${script} \\
        --sample ${sampleID} \\
        --caller clair3_longphase:alignment:${clair3_longphase_vcf} \\
        --caller sniffles:alignment:${sniffles_vcf} \\
        --caller dipcall:assembly:${dipcall_vcf} \\
        --caller hapdiff:assembly:${hapdiff_vcf} \\
        ${trio_args} \\
        --output ${sampleID}.master.vcf
    """
}
