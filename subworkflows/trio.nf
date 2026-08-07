nextflow.enable.dsl=2

include { parse_illumina_reads } from './illumina.nf'
include { atria } from './illumina.nf'
include { bwa_align } from './illumina.nf'
include { filter_bam } from './illumina.nf'
include { deepvariant } from './illumina.nf'
include { YAK_COUNT } from '../modules/yak/main.nf'

// Trio-aware haplotype-family analysis. Children with an entry in the
// family JSON *and* both parents present in the shared Illumina design
// (illumina_design -- the same csv hero_correction reads) get:
//   - a resolved (child, paternal_yak, maternal_yak) row for hifiasm
//     trio-binning at the assembly step
//   - a resolved (child, father, mother, father_vcf, mother_vcf) row for
//     WhatsHap pedigree phasing at the alignment/phasing step
// Parents are never run through the pipeline as full samples -- they only
// need enough to produce a yak db (raw Illumina reads) and a genotype VCF
// (the existing, previously-uncalled atria -> bwa_align -> filter_bam ->
// deepvariant chain). Both downstream consumers (assembly + phasing) reuse
// this single parse/yak/vcf computation rather than repeating it.
workflow trio_family {

    take:
        family_json      // path to JSON: {"child": {"father": "...", "mother": "..."}}, or false/"" to disable
        illumina_design  // shared Illumina design csv (sample,r1,r2)

    main:
        parse_illumina_reads( illumina_design )
        illumina_reads = parse_illumina_reads.out.reads

        def families = family_json ? new groovy.json.JsonSlurper().parse(file(family_json, checkIfExists: true)) : [:]
        trio_rows = families.collect { child, parents -> tuple(child, parents.father, parents.mother) }
        trio_ch = Channel.from(trio_rows)

        parent_ids = (trio_rows.collect { it[1] } + trio_rows.collect { it[2] }).unique()

        parent_reads_ch = illumina_reads.filter { sid, r1, r2 -> sid in parent_ids }

        // paternal/maternal k-mer databases, for hifiasm trio-binning
        YAK_COUNT( parent_reads_ch )
        yak_ch = YAK_COUNT.out.yak_ch

        // genotype-only short-read VCF, for whatshap pedigree phasing
        parent_atria_input = parent_reads_ch.map { sid, r1, r2 -> tuple(sid, [r1, r2]) }
        atria( parent_atria_input )
        bwa_align( atria.out.fqs )
        filter_bam( bwa_align.out.bams )
        deepvariant( filter_bam.out.filter_bam_bai )
        parent_vcf_ch = deepvariant.out.vcfs.map { sid, vcf, tbi -> tuple(sid, vcf) }

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
        trio_assembly = trio_yak_ch      // tuple(child, paternal_yak, maternal_yak)
        trio_phasing  = trio_phasing_ch  // tuple(child, father, mother, father_vcf, mother_vcf)
}
