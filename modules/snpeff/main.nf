process SNPEFF_BUILD {

    tag "${params.snpeff_db}"
    label 'snpeff'

    // published under params.snpeff_data_dir so the next run finds
    // snpEffectPredictor.bin already there and skips rebuilding (see
    // the snpeff_db_exists check in the snpeff subworkflow, subworkflows/ont.nf)
    publishDir "${params.snpeff_data_dir}", mode: params.publish_mode

    input:
        path(genome)
        path(gtf)

    output:
        path("${params.snpeff_db}"), emit: snpeff_db_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    mkdir -p ${params.snpeff_db}
    cp ${genome} ${params.snpeff_db}/sequences.fa
    cp ${gtf} ${params.snpeff_db}/genes.gtf

    snpEff -Xmx${params.snpeff_java_mem} build \\
        -gtf22 \\
        -dataDir \$(pwd) \\
        -configOption ${params.snpeff_db}.genome=${params.snpeff_db} \\
        -v ${params.snpeff_db} \\
        -noCheckCds \\
        -noCheckProtein
    """
}

process SNPEFF {

    tag "$sampleID"
    label 'snpeff'
    errorStrategy 'ignore'

    publishDir "${params.results}/06_variants/snpeff", mode: params.publish_mode

    input:
        tuple val(sampleID), val(vcf)
        path(snpeff_db_dir)

    output:
        tuple val(sampleID), path("${sampleID}.snpeff.vcf"), emit: snpeff_vcf_ch
        tuple val(sampleID), path("${sampleID}.snpeff.csv"), path("${sampleID}.snpeff.html"), optional: true, emit: snpeff_stats_ch

    script:
    // snpEff has no working -t/multi-threading CLI flag (the underlying
    // parallel-stream code exists in SnpEffCmdEff.java but the "-t" case
    // that would enable it is commented out upstream -- true as of the
    // current pcingola/SnpEff master, not just this container's version).
    // csvStats/htmlStats need one coherent pass over the whole VCF, so
    // when stats are requested run the classic single-process path;
    // otherwise get real parallelism by splitting the VCF by chromosome
    // and running one snpEff per chromosome concurrently (bounded to
    // task.cpus), then concatenating -- each chromosome's variants are
    // annotated independently, so this changes nothing about the
    // annotation itself. Each worker still reloads the full genome model
    // regardless of how small its chromosome's slice is, so the per-worker
    // heap (snpeff_scatter_java_mem) is kept modest and the "snpeff" label's
    // --mem is sized for task.cpus concurrent DB loads, not the DB loaded once.
    if (params.snpeff_stats)
        """
        #!/bin/bash
        set -euo pipefail

        snpEff -Xmx${params.snpeff_java_mem} \\
            -dataDir \$(pwd) \\
            -configOption ${params.snpeff_db}.genome=${params.snpeff_db} \\
            -csvStats ${sampleID}.snpeff.csv \\
            -htmlStats ${sampleID}.snpeff.html \\
            ${params.snpeff_args ?: ''} \\
            ${params.snpeff_db} \\
            ${vcf} \\
            > ${sampleID}.snpeff.vcf
        """
    else
        """
        #!/bin/bash
        set -euo pipefail

        grep '^#' ${vcf} > header.vcf

        # single pass: bucket variant lines by chromosome, recording each
        # chromosome's first-appearance order so the concatenated output
        # doesn't need a separate re-sort step
        awk '!/^#/ { if (!seen[\$1]++) print \$1 >> "chroms.txt"; print >> (\$1 ".body.vcf") }' ${vcf}

        > jobs.sh
        while read -r chrom; do
            cat header.vcf "\${chrom}.body.vcf" > "\${chrom}.in.vcf"
            echo "snpEff -Xmx${params.snpeff_scatter_java_mem} -dataDir \$(pwd) -configOption ${params.snpeff_db}.genome=${params.snpeff_db} -noStats ${params.snpeff_args ?: ''} ${params.snpeff_db} \${chrom}.in.vcf > \${chrom}.out.vcf 2> \${chrom}.log" >> jobs.sh
        done < chroms.txt

        xargs -P ${task.cpus} -I{} bash -c '{}' < jobs.sh

        while read -r chrom; do
            [ -s "\${chrom}.out.vcf" ] || { echo "snpEff failed for chromosome \$chrom:" >&2; cat "\${chrom}.log" >&2; exit 1; }
        done < chroms.txt

        first=\$(head -n1 chroms.txt)
        grep '^#' "\${first}.out.vcf" > ${sampleID}.snpeff.vcf
        while read -r chrom; do
            grep -v '^#' "\${chrom}.out.vcf" >> ${sampleID}.snpeff.vcf
        done < chroms.txt
        """
}
