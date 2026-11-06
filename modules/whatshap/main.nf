
process WHATSHAP_PHASE {

    label 'whatshap'
    tag "${sampleID}"
    publishDir "${params.results}/variants/whatshap_phased", mode: params.publish_mode

    input:
        tuple val(ref), val(ref_fai)
        tuple val(sampleID), val(vcf), val(tbi), val(bam), val(bai)

    output:
        tuple val(sampleID), val(bam), val(bai), path("*whatshapphase.vcf.gz"), path("*whatshapphase.vcf.gz.tbi"), emit : whatshap_phase_ch

    script:
    """
    #!/bin/bash

    name=\$(basename ${vcf} .vcf.gz)

    whatshap phase \\
      --reference ${ref} \\
      --output \$name.whatshapphase.vcf.gz \\
      --ignore-read-groups \\
      ${vcf} \\
      ${bam}

    bcftools index -f -t \$name.whatshapphase.vcf.gz

    #mv phased.vcf.gz \$name.whatshapphase.vcf.gz
    #mv phased.vcf.gz.tbi \$name.whatshapphase.vcf.gz.tbi
    """
}

process WHATSHAP_HAPLOTAG {
    
    label 'whatshap'
    tag "${sampleID}"
    publishDir "${params.results}/alignment/whatshap_haplotagged", mode: params.publish_mode

    input:
        tuple val(ref), val(ref_fai)
        tuple val(sampleID), val(bam), val(bai), val(phased_vcf), val(phased_tbi)

    output:
        tuple val(sampleID), path("*.happlotagged.bam"), path("*.happlotagged.bam.bai"), emit : whathap_haplotag_ch

    script:
    """
    #!/bin/bash

    name=\$(basename ${bam} .bam)

    whatshap haplotag \\
      --reference ${ref} \\
      --ignore-read-groups \\
      --output \$name.happlotagged.bam \\
      --output-haplotag-list \$name.happlotagged.tsv \\
      ${phased_vcf} \\
      ${bam}

    samtools index \$name.happlotagged.bam
    """
}

process WHATSHAP_PHASE_TRIO {

    label 'whatshap'
    tag "${child}"
    publishDir "${params.results}/06_variants/whatshap_trio_phased", mode: params.publish_mode

    // Additive, trio-aware refinement of the child's small-variant (SNV +
    // indel) phasing, run alongside (not instead of) the longphase_sv path.
    // Parents contribute genotype-only Mendelian constraints from their own
    // ONT Clair3 VCFs (subworkflows/trio.nf), merged with the child into one
    // multi-sample VCF -- no parental BAM is needed. Sniffles SVs are not
    // included: whatshap only phases SNVs/indels, not symbolic SVs.
    // Outputs: <child>.trio.whatshapphase.vcf.gz (child only, same name as
    // before) and <child>.trio.family.whatshapphase.vcf.gz (all three).
    input:
        tuple val(ref), val(ref_fai)
        tuple val(child), val(father), val(mother), val(child_vcf), val(child_bam), val(child_bai), val(father_vcf), val(mother_vcf)

    output:
        tuple val(child), path("${child}.trio.whatshapphase.vcf.gz"), path("${child}.trio.whatshapphase.vcf.gz.tbi"), emit : whatshap_phase_trio_ch
        tuple val(child), path("${child}.trio.family.whatshapphase.vcf.gz"), path("${child}.trio.family.whatshapphase.vcf.gz.tbi"), emit : whatshap_phase_trio_family_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    printf "%s\\t%s\\t%s\\t%s\\t0\\t0\\n" "${child}_family" "${child}" "${father}" "${mother}" > ${child}.ped
    printf "%s\\t%s\\t0\\t0\\t1\\t0\\n" "${child}_family" "${father}" >> ${child}.ped
    printf "%s\\t%s\\t0\\t0\\t2\\t0\\n" "${child}_family" "${mother}" >> ${child}.ped

    # 1) ONE multi-sample VCF with child + father + mother. In pedigree mode
    #    whatshap only phases the genotypes in its FIRST positional VCF --
    #    every further positional file (VCF or BAM) is read as phasing
    #    EVIDENCE, not as genotypes, so passing the parents' VCFs separately
    #    (as this step used to) gave whatshap no parental genotypes at all and
    #    nothing got phased. Sample columns must match the .ped IDs (CLAIR3
    #    --sample_name). A site absent from a parent's VCF stays missing
    #    (unknown), not assumed hom-ref.
    stage_vcf() {  # stage_vcf <src> <name>: local symlink + index
        ln -sf "\$1" "\$2"
        if [ -e "\$1.tbi" ]; then ln -sf "\$1.tbi" "\$2.tbi"; else bcftools index -f -t "\$2"; fi
    }
    stage_vcf ${child_vcf}  child.vcf.gz
    stage_vcf ${father_vcf} father.vcf.gz
    stage_vcf ${mother_vcf} mother.vcf.gz
    for pair in "child.vcf.gz:${child}" "father.vcf.gz:${father}" "mother.vcf.gz:${mother}"; do
        f=\${pair%%:*}; want=\${pair#*:}
        got=\$(bcftools query -l "\$f" | tr '\\n' ',' | sed 's/,\$//')
        if [ "\$got" != "\$want" ]; then
            echo "[WHATSHAP_PHASE_TRIO] \$f sample column is '\$got', expected '\$want' (family.json ID)" >&2
            exit 1
        fi
    done
    bcftools merge --threads ${task.cpus} -Oz -o family.vcf.gz child.vcf.gz father.vcf.gz mother.vcf.gz
    bcftools index -f -t family.vcf.gz

    # 2) reads must carry a read group whose SM is the child: pedigree mode
    #    assigns reads to samples by @RG SM (--ignore-read-groups is not
    #    allowed with --ped). The alignment step adds no @RG, so add one when
    #    the BAM lacks a matching SM -- otherwise no reads count for the child.
    #    pysam (always present -- whatshap depends on it) runs the samtools
    #    commands, so this doesn't rely on a samtools binary in the image.
    if python3 -c "
import pysam, sys
rgs = pysam.AlignmentFile('${child_bam}').header.to_dict().get('RG', [])
sys.exit(0 if any(rg.get('SM') == '${child}' for rg in rgs) else 1)
"; then
        reads_bam=${child_bam}
    else
        echo "[WHATSHAP_PHASE_TRIO] ${child_bam} has no @RG SM:${child} -- adding one" >&2
        python3 -c "
import pysam
pysam.addreplacerg('-@', '${task.cpus}', '-r', 'ID:${child}', '-r', 'SM:${child}', '-m', 'overwrite_all',
                   '-o', '${child}.rg.bam', '${child_bam}', catch_stdout=False)
pysam.index('-@', '${task.cpus}', '${child}.rg.bam')
"
        reads_bam=${child}.rg.bam
    fi

    whatshap phase \\
      --ped ${child}.ped \\
      --reference ${ref} \\
      --output ${child}.trio.family.whatshapphase.vcf.gz \\
      family.vcf.gz \\
      \$reads_bam

    bcftools index -f -t ${child}.trio.family.whatshapphase.vcf.gz
    bcftools view -s ${child} -Oz -o ${child}.trio.whatshapphase.vcf.gz ${child}.trio.family.whatshapphase.vcf.gz
    bcftools index -f -t ${child}.trio.whatshapphase.vcf.gz
    rm -f ${child}.rg.bam ${child}.rg.bam.bai

    # 3) say how much actually got phased -- a silent all-unphased output is
    #    exactly how the previous version failed
    bcftools query -i 'TYPE="snp" && GT="het"' -f '[%GT]\\n' ${child}.trio.whatshapphase.vcf.gz \\
        | awk '{ if (\$1 ~ /\\|/) p++; n++ } END { printf "[WHATSHAP_PHASE_TRIO] ${child}: %d / %d het SNVs phased (%.1f%%)\\n", p, n, n ? 100*p/n : 0 }' >&2
    """
}