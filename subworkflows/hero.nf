nextflow.enable.dsl=2

include { SORT_READS_FOR_ROPEBWT2 } from '../modules/ropebwt2/main.nf'
include { ROPEBWT2 } from '../modules/ropebwt2/main.nf'
include { FMLRC2_CONVERT } from '../modules/fmlrc2/main.nf'
include { FMLRC2_CORRECT } from '../modules/fmlrc2/main.nf'
include { HERO } from '../modules/hero/main.nf'
include { SAMTOOLS_FASTQ_TO_BAM } from '../modules/samtools/main.nf'
include { parse_illumina_reads } from './illumina.nf'
include { atria } from './illumina.nf'

// Hybrid error correction for samples that have matched Illumina reads;
// samples without a matching entry in the Illumina design pass through
// unchanged. Runs on the post-dorado-trim fastq (not right after
// parse_design) because DORADO_TRIM only ever reads the unaligned bam --
// the merged fastq from parse_design is never consumed downstream, so
// correcting it there would be a dead end.
workflow hero_correction {

    take:
        ont_ch           // tuple(sampleID, unal_bam, fastq)
        illumina_design  // path to illumina design csv (sample,r1,r2), or false/empty to disable

    main:
        parse_illumina_reads( illumina_design )

        // adapter-trim unconditionally -- input fastqs aren't guaranteed to
        // already be trimmed, and re-trimming already-clean reads is a no-op
        atria( parse_illumina_reads.out.reads.map { sampleID, r1, r2 -> tuple(sampleID, [r1, r2]) } )
        illumina_reads = atria.out.fqs.map { sampleID, fqs -> tuple(sampleID, fqs[0], fqs[1]) }

        gated = ont_ch
            .join(illumina_reads, remainder: true)
            .branch {
                has_illumina: it[3] != null
                no_illumina: true
            }

        to_correct = gated.has_illumina
        skip_correction = gated.no_illumina
            // join(remainder:true) pads an unmatched item with a single null,
            // not one null per missing right-hand field -- so this tuple is
            // (sampleID, bam, fastq, null), not the 5-element shape matched
            // items have. Index instead of destructuring by name so it works
            // for both.
            .map { row -> tuple(row[0], row[1], row[2]) }

        SORT_READS_FOR_ROPEBWT2( to_correct.map { sampleID, bam, fastq, r1, r2 -> tuple(sampleID, r1, r2) } )
        ROPEBWT2( SORT_READS_FOR_ROPEBWT2.out.sorted_ch )

        FMLRC2_CONVERT( ROPEBWT2.out.ropebwt2_ch )

        fmlrc2_correct_input = FMLRC2_CONVERT.out.fmlrc2_convert_ch
            .join( to_correct.map { sampleID, bam, fastq, r1, r2 -> tuple(sampleID, fastq) } )

        FMLRC2_CORRECT( fmlrc2_correct_input )

        hero_input = FMLRC2_CORRECT.out.fmlrc2_correct_ch
            .join( to_correct.map { sampleID, bam, fastq, r1, r2 -> tuple(sampleID, r1, r2) } )

        HERO( hero_input )

        rebuild_input = HERO.out.hero_ch
            .join( to_correct.map { sampleID, bam, fastq, r1, r2 -> tuple(sampleID, bam) } )

        SAMTOOLS_FASTQ_TO_BAM( rebuild_input )

        ont_ch_final = SAMTOOLS_FASTQ_TO_BAM.out.fastq_to_bam_ch
            .mix( skip_correction )

    emit:
        ont_ch = ont_ch_final
}
