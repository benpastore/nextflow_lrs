nextflow.enable.dsl=2

// Mild vs severe haplotype comparison -> genetic models. See modules/haplotype/main.nf.
include { SIB_HAPLOTYPE_SHARING; HAPLOTYPE_MODIFIER } from '../modules/haplotype/main.nf'

workflow haplotype_analysis {

    take :
        family_json
        parent_vcf       // tuple(parentID, clair3 vcf)       -- trio_family.out.parent_vcf
        child_vcf        // tuple(childID, phased vcf, tbi)    -- WHATSHAP_PHASE_TRIO's child VCF
        vcf_manifest     // childID<TAB>snpEff VCF name (consolidated, with INFO/ORIGIN)
        vcfs             // collected snpEff VCFs
        kg_af_vcf        // or file("NO_FILE")

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
            kg_af_vcf,
            file("${params.bin}/_modifier_common.py"),
            file("${params.bin}/haplotype_modifier.py")
        )

    emit :
        recessive_modifier_ch = HAPLOTYPE_MODIFIER.out.recessive_modifier_ch
        dominant_modifier_ch  = HAPLOTYPE_MODIFIER.out.dominant_modifier_ch
        haplotype_variants_ch = HAPLOTYPE_MODIFIER.out.haplotype_variants_ch
}
