
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


process WHATSHAP_PHASE_FAMILY {

    label 'whatshap'
    tag "${family}"
    publishDir "${params.results}/06_variants/whatshap_family_phased", mode: params.publish_mode

    // Joint pedigree phasing of a whole family -- every child + both parents
    // in ONE whatshap run (one .ped, one multi-sample VCF, every member's
    // reads), then split back into one phased VCF per child. With
    // --run_trio_analysis the per-child output REPLACES longphase_sv as the
    // child's phased small-variant VCF for consolidation and the haplotype
    // analysis (see main.nf stage 5).
    //
    // Compared with phasing each child as its own trio, the parents are
    // themselves phased (their reads + transmission to every child), every
    // child is phased against the same parental haplotypes, and a genotype
    // that's Mendelian-inconsistent in one child is judged against all of
    // them. Parents' reads are optional (params.whatshap_family_parent_reads);
    // without them the parents contribute genotypes only.
    //
    // Genotypes: children's longphase_sv VCFs + parents' Clair3 VCFs, merged.
    // Sniffles SVs aren't phased (whatshap phases SNVs/indels only).
    // Run per chromosome in parallel inside this job
    // (params.whatshap_family_jobs at a time), then concatenated.
    //
    // Outputs: <child>.trio.whatshapphase.vcf.gz per child (child column,
    // sites where the child carries ALT -- same name/content contract as the
    // old per-child trio step) and <family>.family.whatshapphase.vcf.gz.
    input:
        tuple val(ref), val(ref_fai)
        tuple val(family), val(father), val(mother),
            val(father_vcf), val(mother_vcf),
            val(father_bam), val(mother_bam),
            val(children), val(child_vcfs), val(child_bams)

    output:
        tuple path("*.trio.whatshapphase.vcf.gz"), path("*.trio.whatshapphase.vcf.gz.tbi"), emit : child_vcfs_ch
        tuple val(family), path("${family}.family.whatshapphase.vcf.gz"), path("${family}.family.whatshapphase.vcf.gz.tbi"), emit : family_vcf_ch
        tuple val(family), path("${family}.recombination.tsv"), emit : recombination_ch

    script:
    def kids     = children instanceof List ? children : [children]
    def kid_vcfs = child_vcfs instanceof List ? child_vcfs : [child_vcfs]
    def kid_bams = child_bams instanceof List ? child_bams : [child_bams]
    def ped_rows = kids.collect { c -> "${family}\\t${c}\\t${father}\\t${mother}\\t0\\t0" }.join("\\n")
    // "sample:vcf" / "sample:bam" pairs, children first
    def vcf_pairs = ([kids, kid_vcfs].transpose().collect { c, v -> "${c}:${v}" } +
                     ["${father}:${father_vcf}", "${mother}:${mother_vcf}"]).join(' ')
    def bam_pairs = ([kids, kid_bams].transpose().collect { c, b -> "${c}:${b}" } +
                     (father_bam != 'NONE' ? ["${father}:${father_bam}"] : []) +
                     (mother_bam != 'NONE' ? ["${mother}:${mother_bam}"] : [])).join(' ')
    def n_jobs = Math.max(1, Math.min(task.cpus as int, (params.whatshap_family_jobs ?: 6) as int))
    """
    #!/bin/bash
    set -euo pipefail

    # 1) pedigree: every child of this father/mother pair
    printf "${ped_rows}\\n" > ${family}.ped
    printf "%s\\t%s\\t0\\t0\\t1\\t0\\n" "${family}" "${father}" >> ${family}.ped
    printf "%s\\t%s\\t0\\t0\\t2\\t0\\n" "${family}" "${mother}" >> ${family}.ped

    # 2) ONE multi-sample VCF (whatshap's pedigree mode phases only the
    #    genotypes in its first positional VCF; further positional files are
    #    read as phasing evidence). Sample columns must match the .ped IDs.
    #    A site absent from a member's VCF stays missing (unknown).
    merge_list=()
    for pair in ${vcf_pairs}; do
        s=\${pair%%:*}; f=\${pair#*:}
        ln -sf "\$f" "gt_\$s.vcf.gz"
        if [ -e "\$f.tbi" ]; then ln -sf "\$f.tbi" "gt_\$s.vcf.gz.tbi"; else bcftools index -f -t "gt_\$s.vcf.gz"; fi
        got=\$(bcftools query -l "gt_\$s.vcf.gz" | tr '\\n' ',' | sed 's/,\$//')
        if [ "\$got" != "\$s" ]; then
            echo "[WHATSHAP_PHASE_FAMILY] \$f sample column is '\$got', expected '\$s' (family.json ID)" >&2
            exit 1
        fi
        merge_list+=("gt_\$s.vcf.gz")
    done
    bcftools merge --threads ${task.cpus} -Oz -o family.vcf.gz "\${merge_list[@]}"
    bcftools index -f -t family.vcf.gz

    # 3) reads: pedigree mode assigns reads to samples by @RG SM, so every
    #    BAM needs a read group whose SM is its sample; add one when missing
    bams=()
    for pair in ${bam_pairs}; do
        s=\${pair%%:*}; b=\${pair#*:}
        if python3 -c "
import pysam, sys
rgs = pysam.AlignmentFile('\$b').header.to_dict().get('RG', [])
sys.exit(0 if any(rg.get('SM') == '\$s' for rg in rgs) else 1)
"; then
            bams+=("\$b")
        else
            echo "[WHATSHAP_PHASE_FAMILY] \$b has no @RG SM:\$s -- adding one" >&2
            python3 -c "
import pysam
pysam.addreplacerg('-@', '${task.cpus}', '-r', 'ID:\$s', '-r', 'SM:\$s', '-m', 'overwrite_all',
                   '-o', 'rg_\$s.bam', '\$b', catch_stdout=False)
pysam.index('-@', '${task.cpus}', 'rg_\$s.bam')
"
            bams+=("rg_\$s.bam")
        fi
    done

    # 4) phase per chromosome, ${n_jobs} at a time. Each chromosome's input
    #    is subset first so its output holds only that chromosome.
    mkdir -p chunks
    phase_chrom() {
        c=\$1
        bcftools view -r "\$c" -Oz -o "chunks/in.\$c.vcf.gz" family.vcf.gz
        bcftools index -f -t "chunks/in.\$c.vcf.gz"
        # --recombination-list crashes in whatshap 2.1 (AssertionError in
        # find_recombination) on a chromosome with nothing phased, e.g. chrM;
        # retry such a chromosome without it (no recombination rows for it)
        if whatshap phase \\
              --ped ${family}.ped \\
              --reference ${ref} \\
              --chromosome "\$c" \\
              --output "chunks/out.\$c.vcf.gz" \\
              --recombination-list "chunks/recomb.\$c.tsv" \\
              "chunks/in.\$c.vcf.gz" \\
              "\${bams[@]}" > "chunks/\$c.log" 2>&1; then
            touch "chunks/\$c.done"
        else
            echo "[WHATSHAP_PHASE_FAMILY] \$c failed with --recombination-list; retrying without it" >&2
            rm -f "chunks/out.\$c.vcf.gz" "chunks/recomb.\$c.tsv"
            whatshap phase \\
              --ped ${family}.ped \\
              --reference ${ref} \\
              --chromosome "\$c" \\
              --output "chunks/out.\$c.vcf.gz" \\
              "chunks/in.\$c.vcf.gz" \\
              "\${bams[@]}" > "chunks/\$c.log" 2>&1 \\
              && touch "chunks/\$c.done"
        fi
    }
    chroms=\$(bcftools index -s family.vcf.gz | cut -f1)
    for c in \$chroms; do
        while [ "\$(jobs -rp | wc -l)" -ge ${n_jobs} ]; do sleep 5; done
        phase_chrom "\$c" &
    done
    wait || true
    failed=0
    for c in \$chroms; do
        if [ ! -e "chunks/\$c.done" ]; then
            echo "[WHATSHAP_PHASE_FAMILY] whatshap failed on \$c:" >&2
            tail -20 "chunks/\$c.log" >&2 || true
            failed=1
        fi
    done
    [ \$failed -eq 0 ] || exit 1

    bcftools concat -Oz -o ${family}.family.whatshapphase.vcf.gz \$(for c in \$chroms; do echo "chunks/out.\$c.vcf.gz"; done)
    bcftools index -f -t ${family}.family.whatshapphase.vcf.gz

    # 4b) whatshap's inferred crossovers (--recombination-list), one table
    #     for the family: the first chunk's header line(s), then every row.
    #     Exploratory -- nothing downstream reads it yet.
    : > ${family}.recombination.tsv
    hdr=0
    for c in \$chroms; do
        f="chunks/recomb.\$c.tsv"
        [ -s "\$f" ] || continue
        if [ \$hdr -eq 0 ]; then grep '^#' "\$f" >> ${family}.recombination.tsv || true; hdr=1; fi
        grep -v '^#' "\$f" >> ${family}.recombination.tsv || true
    done
    echo "[WHATSHAP_PHASE_FAMILY] \$(grep -vc '^#' ${family}.recombination.tsv || true) recombination events in ${family}.recombination.tsv" >&2

    # 5) split: child column only, sites where the child carries ALT
    for child in ${kids.join(' ')}; do
        bcftools view -s "\$child" ${family}.family.whatshapphase.vcf.gz \\
            | bcftools view -i 'GT="alt"' -Oz -o "\$child.trio.whatshapphase.vcf.gz"
        bcftools index -f -t "\$child.trio.whatshapphase.vcf.gz"
    done

    # 6) how much got phased, per member
    for s in ${kids.join(' ')} ${father} ${mother}; do
        bcftools view -s "\$s" ${family}.family.whatshapphase.vcf.gz \\
            | bcftools query -i 'TYPE="snp" && GT="het"' -f '[%GT]\\n' \\
            | awk -v s="\$s" '{ if (\$1 ~ /\\|/) p++; n++ } END { printf "[WHATSHAP_PHASE_FAMILY] %s: %d / %d het SNVs phased (%.1f%%)\\n", s, p, n, n ? 100*p/n : 0 }' >&2
    done
    rm -rf chunks rg_*.bam rg_*.bam.bai
    """
}
