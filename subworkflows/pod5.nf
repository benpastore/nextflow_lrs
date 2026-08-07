nextflow.enable.dsl=2

include { DORADO_BASECALL } from '../modules/dorado/main.nf'
include { SAMTOOLS_CONVERT_BAM_TO_FASTQ } from '../modules/samtools/main.nf'

// Alternative pipeline entry point: basecall directly from pod5 instead of
// starting from an already-basecalled unaligned bam/fastq (parse_design).
// Produces the same (sampleID, fastq, unal_bam) shape parse_design.out.ont
// does, so it drops straight into DORADO_TRIM in main.nf unchanged.
workflow pod5_basecall {

    take:
        pod5_design  // path to csv (sample,pod5) -- pod5 column is a directory of pod5 files per sample

    main:
        pod5_reads = Channel
            .fromPath(file(pod5_design, checkIfExists: true))
            .splitCsv(header: ['sample', 'pod5'], sep: ',', skip: 1)
            .map { row -> tuple(row.sample, file(row.pod5)) }

        DORADO_BASECALL( pod5_reads )

        SAMTOOLS_CONVERT_BAM_TO_FASTQ( DORADO_BASECALL.out.dorado_basecall_ch )

        ont_ch = SAMTOOLS_CONVERT_BAM_TO_FASTQ.out.samtool_convert_unalbam2fastq
            .map { sampleID, bam, fastq -> tuple(sampleID, fastq, bam) }

    emit:
        ont = ont_ch  // tuple(sampleID, fastq, unal_bam)
}
