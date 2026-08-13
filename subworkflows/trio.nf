nextflow.enable.dsl=2

include { YAK_COUNT } from '../modules/yak/main.nf'
include { minimap2 } from './ont.nf'
include { clair3 } from './ont.nf'

// Trio-aware haplotype-family analysis. Parents have no Illumina data --
// they're sequenced on ONT the same as children, just listed as ordinary
// rows in the same design.csv/pod5_design. Since they're never meant to be
// run through the pipeline as full samples (no assembly, SV calling,
// phasing, paraphase, or modkit), this subworkflow splits them out of the
// shared ont_ch as soon as they're identifiable (family_json), gives them
// just enough processing to be useful for trio analysis --
//   - yak k-mer database (straight from their own ONT reads) for hifiasm
//     trio-binning
//   - minimap2 + clair3 genotype-only VCF (the pipeline's own ONT
//     small-variant caller, reused as-is) for WhatsHap pedigree phasing
// -- and hands back non_parent_ont_ch for main.nf to use in place of the
// full ont_ch for every other section of the pipeline.
workflow trio_family {

    take:
        family_json  // path to JSON: {"child": {"father": "...", "mother": "..."}}, or false/"" to disable
        ont_ch        // tuple(sampleID, unal_bam, fastq) -- every sample (children + parents), post-HERO
        reference_index    // minimap2 .mmi index
        samtools_fai_index // tuple(ref, ref_fai)

    main:
        def families = family_json ? new groovy.json.JsonSlurper().parse(file(family_json, checkIfExists: true)) : [:]
        trio_rows = families.collect { child, parents -> tuple(child, parents.father, parents.mother) }
        trio_ch = Channel.from(trio_rows)

        parent_ids = (trio_rows.collect { it[1] } + trio_rows.collect { it[2] }).unique()

        parent_ont_ch = ont_ch.filter { sampleID, bam, fastq -> sampleID in parent_ids }
        non_parent_ont_ch = ont_ch.filter { sampleID, bam, fastq -> !(sampleID in parent_ids) }

        // paternal/maternal k-mer databases, for hifiasm trio-binning
        YAK_COUNT( parent_ont_ch.map { sampleID, bam, fastq -> tuple(sampleID, [fastq]) } )
        yak_ch = YAK_COUNT.out.yak_ch

        // genotype-only VCF, for whatshap pedigree phasing -- same
        // minimap2/clair3 the rest of the pipeline uses, just run on
        // parents only and not carried any further (no sniffles/longphase/
        // straglr/spectre/paraphase/modkit for parents)
        minimap2( reference_index, parent_ont_ch.map { sampleID, bam, fastq -> tuple(sampleID, fastq) } )
        clair3( samtools_fai_index, minimap2.out.bams )
        parent_vcf_ch = clair3.out.clair3_ch.map { sampleID, vcf, tbi -> tuple(sampleID, vcf) }

        // resolve full trios against yak (join keyed on father, then remapped and
        // joined again keyed on mother -- Nextflow's join() only keys on one field
        // at a time, so a 2-parent join needs this double-remap-and-join)
        trio_yak_ch = trio_ch
            .map { child, father, mother -> tuple(father, child, mother) }
            .join( yak_ch )
            .map { father, child, mother, pat_yak -> tuple(mother, child, father, pat_yak) }
            .join( yak_ch )
            .map { mother, child, father, pat_yak, mat_yak -> tuple(child, pat_yak, mat_yak) }

        // resolve full trios against parent vcfs, same double-join pattern
        trio_phasing_ch = trio_ch
            .map { child, father, mother -> tuple(father, child, mother) }
            .join( parent_vcf_ch )
            .map { father, child, mother, father_vcf -> tuple(mother, child, father, father_vcf) }
            .join( parent_vcf_ch )
            .map { mother, child, father, father_vcf, mother_vcf -> tuple(child, father, mother, father_vcf, mother_vcf) }

    emit:
        trio_assembly     = trio_yak_ch        // tuple(child, paternal_yak, maternal_yak)
        trio_phasing       = trio_phasing_ch    // tuple(child, father, mother, father_vcf, mother_vcf)
        non_parent_ont_ch = non_parent_ont_ch  // tuple(sampleID, unal_bam, fastq) -- everyone except parents
}
