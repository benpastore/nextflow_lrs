nextflow.enable.dsl=2

include { HERRO_PREPROCESS } from '../modules/herro/main.nf'
include { HERRO_ALIGN_BATCHES } from '../modules/herro/main.nf'
include { HERRO_INFERENCE } from '../modules/herro/main.nf'
include { SAMTOOLS_FASTQ_TO_BAM } from '../modules/samtools/main.nf'

// GPU deep-learning ONT self-correction (lbcb-sci/herro), replacing the old
// hybrid Illumina+ONT "HERO" tool. herro has no notion of matched Illumina
// data, so unlike the old illumina_design-gated path, every sample in
// ont_ch goes through the same three steps uniformly -- no branch/skip.
workflow herro_correction {

    take:
        ont_ch  // tuple(sampleID, unal_bam, fastq) -- every sample, post-dorado-trim

    main:
        HERRO_PREPROCESS( ont_ch )

        HERRO_ALIGN_BATCHES( HERRO_PREPROCESS.out.preprocessed_ch )

        HERRO_INFERENCE( HERRO_ALIGN_BATCHES.out.batches_ch )

        // reimport corrected fastq into an unaligned bam, reusing the
        // original bam's @RG (needed later for dorado polish's model
        // auto-resolution -- see SAMTOOLS_FASTQ_TO_BAM)
        rebuild_input = HERRO_INFERENCE.out.herro_ch
            .join( ont_ch.map { sampleID, bam, fastq -> tuple(sampleID, bam) } )

        SAMTOOLS_FASTQ_TO_BAM( rebuild_input )

    emit:
        ont_ch = SAMTOOLS_FASTQ_TO_BAM.out.fastq_to_bam_ch
}
