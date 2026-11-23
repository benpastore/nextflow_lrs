nextflow.enable.dsl=2

// Mild vs severe sibling haplotype sharing with hap-ibd. See modules/haplotype/main.nf.
include { ORIENT_PHASE; HAP_IBD; IBD_SEGMENTS; MILD_HAPLOTYPE_VARIANTS } from '../modules/haplotype/main.nf'

workflow haplotype_analysis {

    take :
        family_json
        family_vcf       // tuple(family, family.whatshapphase.vcf.gz, tbi) -- WHATSHAP_PHASE_FAMILY
        samtools_fai_index

    main :
        // 1. orient each child's phase paternal|maternal, children only
        ORIENT_PHASE( family_vcf, family_json, file("${params.bin}/orient_family_phase.py") )

        // 2. hap-ibd per family (haplotype 1 = paternal, 2 = maternal)
        HAP_IBD( ORIENT_PHASE.out.oriented_ch, samtools_fai_index )

        // 3. per mild x severe pair, per parent: shared / not_shared, IBD0/1/2
        fam_ch = ORIENT_PHASE.out.oriented_ch.join( HAP_IBD.out.ibd_ch )
        IBD_SEGMENTS(
            family_json,
            fam_ch.map { fam, vcf, ibd -> fam }.collect(),
            fam_ch.map { fam, vcf, ibd -> vcf }.collect(),
            fam_ch.map { fam, vcf, ibd -> ibd }.collect(),
            file("${params.bin}/hapibd_segments.py")
        )

    emit :
        parent_segments_ch = IBD_SEGMENTS.out.parent_segments_ch   // per parent: shared / not_shared
        ibd_segments_ch    = IBD_SEGMENTS.out.ibd_segments_ch      // IBD0/1/2
        ibd_summary_ch     = IBD_SEGMENTS.out.ibd_summary_ch       // genome fractions per pair
}

// Mild-child variants on a parental copy the severe sibling didn't inherit.
workflow mild_haplotype_variants {

    take :
        family_json
        variant_vcf      // tuple(childID, annotated master VCF)
        parent_segments  // cohort.sib_parent_segments.tsv

    main :
        manifest = variant_vcf
            .map { sid, vcf -> "${sid}\t${vcf.name}" }
            .collectFile(name: 'mild_variants_vcf_manifest.tsv', newLine: true, sort: true)
        MILD_HAPLOTYPE_VARIANTS(
            family_json, manifest, variant_vcf.map { sid, vcf -> vcf }.collect(),
            parent_segments, file("${params.bin}/mild_haplotype_variants.py")
        )

    emit :
        variants_ch = MILD_HAPLOTYPE_VARIANTS.out.variants_ch
        genes_ch    = MILD_HAPLOTYPE_VARIANTS.out.genes_ch
}

// Variant-based modifier search (one cohort task): per discordant family, the
// mild child's passing genes (homozygous / compound_het / regulatory) vs its
// severe sib(s), concordant severe children as extra exclusion controls, then the
// cohort count of mild children per gene; for the recurrent genes every mild
// and severe variant in the gene goes to AlphaGenome / AlphaMissense as two
// deduplicated VCFs. See modules/modifier/main.nf.
include { MODIFIER_VARIANTS; GENE_RECURRENCE } from '../modules/modifier/main.nf'

workflow modifier_analysis {

    take :
        family_json      // path to family.json
        variant_vcf      // tuple(childID, annotated master VCF: ORIGIN + ANN + AF)
        ibd_segments     // cohort.sib_ibd_segments.tsv (hapibd_segments.py)
        gtf              // gene model on the reference (protein-coding = has CDS)

    main :
        // one task for the cohort: every mild / severe child goes in;
        // modifier_variants.py groups them into families (father + mother),
        // runs the rules per discordant family and uses the severe children
        // of concordant (severe-only) families as extra exclusion controls
        def fam = new groovy.json.JsonSlurper().parseText(family_json.text)
        child_ch = variant_vcf
            .filter { sid, vcf -> fam[sid] && fam[sid].phenotype in ['mild', 'severe'] }
            .toSortedList { a, b -> a[0] <=> b[0] }
        sample_args = child_ch.map { kids -> kids.collect { sid, vcf -> "-sample_vcf ${sid}=${vcf.name}" }.join(' ') }
        vcf_files = child_ch.map { kids -> kids.collect { sid, vcf -> vcf } }

        MODIFIER_VARIANTS(
            sample_args, vcf_files, family_json, ibd_segments, gtf,
            file("${params.bin}/_modifier_common.py"),
            file("${params.bin}/_cnv_repeat_common.py"),
            file("${params.bin}/modifier_variants.py")
        )

        GENE_RECURRENCE(
            MODIFIER_VARIANTS.out.genes_ch.collect(),
            MODIFIER_VARIANTS.out.gene_variants_ch.collect(),
            file("${params.bin}/gene_recurrence.py")
        )

    emit :
        passing_genes_ch        = GENE_RECURRENCE.out.passing_genes_ch
        recurrence_ch           = GENE_RECURRENCE.out.recurrence_ch
        // AlphaGenome / AlphaMissense input: every mild / severe variant in the recurrent genes
        mild_candidate_vcf_ch   = GENE_RECURRENCE.out.mild_candidate_vcf_ch
        severe_candidate_vcf_ch = GENE_RECURRENCE.out.severe_candidate_vcf_ch
}

// Mild vs severe promoter methylation per parental haplotype. See
// modules/methylation_compare/main.nf.
include { METHYLATION_COMPARE } from '../modules/methylation_compare/main.nf'

workflow methylation_compare {

    take :
        family_json
        modkit_beds      // tuple(childID, hp1.bed.gz, hp2.bed.gz, combined.bed.gz)
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
