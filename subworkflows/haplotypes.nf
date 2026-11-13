nextflow.enable.dsl=2

// Mild vs severe haplotype comparison -> genetic models. See modules/haplotype/main.nf.
include { SIB_HAPLOTYPE_SHARING; HAPLOTYPE_MODIFIER } from '../modules/haplotype/main.nf'

workflow haplotype_analysis {

    take :
        family_json
        parent_vcf       // tuple(parentID, clair3 vcf)       -- trio_family.out.parent_vcf
        child_vcf        // tuple(childID, phased vcf, tbi)    -- WHATSHAP_PHASE_FAMILY child VCF
        vcf_manifest     // childID<TAB>snpEff VCF name (consolidated, with INFO/ORIGIN)
        vcfs             // collected snpEff VCFs (with INFO/gnomAD_AF_NFE, KG_AF_EUR)

    main :
        // 1. which parental haplotypes each mild/severe pair shares
        geno_vcf = parent_vcf.mix( child_vcf.map { sid, vcf, tbi -> tuple(sid, vcf) } )
        geno_manifest = geno_vcf
            .map { sid, vcf -> "${sid}\t${vcf.name}" }
            .collectFile(name: 'haplotype_vcf_manifest.tsv', newLine: true, sort: true)
        SIB_HAPLOTYPE_SHARING(
            family_json, geno_manifest, geno_vcf.map { sid, vcf -> vcf }.collect(),
            file("${params.bin}/sib_haplotype_sharing.py")
        )

        // 2. mild variants on shared vs mild-specific haplotypes -> models
        HAPLOTYPE_MODIFIER(
            family_json, vcf_manifest, vcfs,
            SIB_HAPLOTYPE_SHARING.out.parent_segments_ch,
            file("${params.bin}/_modifier_common.py"),
            file("${params.bin}/haplotype_modifier.py")
        )

    emit :
        recessive_modifier_ch = HAPLOTYPE_MODIFIER.out.recessive_modifier_ch
        dominant_modifier_ch  = HAPLOTYPE_MODIFIER.out.dominant_modifier_ch
        haplotype_variants_ch = HAPLOTYPE_MODIFIER.out.haplotype_variants_ch
        parent_segments_ch    = SIB_HAPLOTYPE_SHARING.out.parent_segments_ch
}

// Mild vs severe promoter methylation per parental haplotype. See
// modules/methylation_compare/main.nf.
include { METHYLATION_COMPARE } from '../modules/methylation_compare/main.nf'

workflow methylation_compare {

    take :
        family_json
        modkit_beds      // tuple(childID, hp1.bed.gz, hp2.bed.gz, ungrouped.bed.gz)
        phased_vcf       // tuple(childID, longphase_sv vcf, tbi) -- source of the HP tags
        parent_vcf       // tuple(parentID, clair3 vcf)
        parent_segments  // cohort.sib_parent_segments.tsv
        gtf

    main :
        child_ch = modkit_beds.join( phased_vcf.map { sid, vcf, tbi -> tuple(sid, vcf) } )
        manifest = child_ch
            .map { sid, hp1, hp2, ung, vcf -> "${sid}\t${hp1.name}\t${hp2.name}\t${ung.name}\t${vcf.name}" }
            .collectFile(name: 'methylation_manifest.tsv', newLine: true, sort: true)
        child_files = child_ch.flatMap { sid, hp1, hp2, ung, vcf -> [hp1, hp2, ung, vcf] }.collect()

        parent_manifest = parent_vcf
            .map { sid, vcf -> "${sid}\t${vcf.name}" }
            .collectFile(name: 'methylation_parent_manifest.tsv', newLine: true, sort: true)

        METHYLATION_COMPARE(
            family_json, manifest, child_files,
            parent_manifest, parent_vcf.map { sid, vcf -> vcf }.collect(),
            parent_segments, gtf,
            file("${params.bin}/_modifier_common.py"),
            file("${params.bin}/_cnv_repeat_common.py"),
            file("${params.bin}/methylation_compare.py")
        )

    emit :
        promoter_methylation_ch = METHYLATION_COMPARE.out.promoter_methylation_ch
}
