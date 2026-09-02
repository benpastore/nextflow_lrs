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

        // resolve full trios against yak/parent-vcf by collecting each into a
        // plain map and looking children up against it, rather than
        // Channel.join() -- join() consumes each channel item once, so it
        // silently drops any child past the first one that shares a parent
        // (e.g. two children of the same father/mother both keying on that
        // father: only the first join match gets paired, the second child
        // just vanishes from the channel with no error). Parent counts are
        // always small, so collecting to a map is cheap and correct for
        // however many children reference the same parent.
        yak_map_ch = yak_ch
            .toList()
            .map { list -> list.collectEntries { sampleID, yak -> [(sampleID): yak] } }

        trio_yak_ch = trio_ch
            .combine( yak_map_ch )
            .map { child, father, mother, yak_map -> tuple(child, yak_map[father], yak_map[mother]) }

        vcf_map_ch = parent_vcf_ch
            .toList()
            .map { list -> list.collectEntries { sampleID, vcf -> [(sampleID): vcf] } }

        trio_phasing_ch = trio_ch
            .combine( vcf_map_ch )
            .map { child, father, mother, vcf_map -> tuple(child, father, mother, vcf_map[father], vcf_map[mother]) }

    emit:
        trio_assembly     = trio_yak_ch        // tuple(child, paternal_yak, maternal_yak)
        trio_phasing      = trio_phasing_ch    // tuple(child, father, mother, father_vcf, mother_vcf)
        non_parent_ont_ch = non_parent_ont_ch  // tuple(sampleID, unal_bam, fastq) -- everyone except parents
}
