#!/usr/bin/env nextflow

include { ATRIA } from '../modules/atria/main.nf'
include { BWA_INDEX } from '../modules/bwa/main.nf'
include { BWA_MEM } from '../modules/bwa/main.nf'
include { FILTER_BAM } from '../modules/samtools/main.nf'
include { INDEX_BAM } from '../modules/samtools/main.nf'
include { BAM_TO_BW } from '../modules/deeptools/main.nf'
include { DEEPVARIANT_CALL_VARIANTS } from '../modules/deepvariant/main.nf'

nextflow.enable.dsl=2

// Shared parser for the "sample,r1,r2" Illumina design csv -- used both by
// hero_correction (children needing hybrid error correction) and
// trio_family (parents needing a yak db / short-read VCF), so any sampleID
// with Illumina reads only has to be declared once, in one file.
workflow parse_illumina_reads {

    take:
        illumina_design  // path to csv (sample,r1,r2), or false/"" to disable

    main:
        if (illumina_design) {
            reads_ch = Channel
                .fromPath(file(illumina_design, checkIfExists: true))
                .splitCsv(header: ['sample', 'r1', 'r2'], sep: ',', skip: 1)
                .map { row -> tuple(row.sample, file(row.r1), file(row.r2)) }
        } else {
            reads_ch = Channel.empty()
        }

    emit:
        reads = reads_ch
}

workflow atria {

    take : 
        data

    main :
        /*
        * Trim reads of adapters and low quality sequences
        */
        ATRIA( data )

        fq_ch = ATRIA.out.fq_ch

        if ( params.fastqc ){
            FASTQC( data )
        }
    
    emit : 
        fqs = fq_ch
}

workflow bwa_align {

    take : data

    main :
        genome_fasta = file("${params.genome}")
        genome_name = "${genome_fasta.baseName}"
        bwa_index_path = "${params.index}/bwa/${genome_name}"
        bwa_index = "${bwa_index_path}/${genome_name}"
        bwa_chrom_sizes = "${bwa_index_path}/chrom_sizes.txt"

        bwa_exists = file(bwa_index_path).exists()

        if (bwa_exists == true){
            bwa_build = false
        } else {
            bwa_build = true
        }

        if (bwa_build == true){
            BWA_INDEX( params.genome, bwa_index )
            bwa_index_ch = BWA_INDEX.out.bwa_index_ch
        } else {
            bwa_index_ch = bwa_index 
        }

        BWA_MEM( bwa_index_ch, data )
    
    emit :
        bams = BWA_MEM.out.bwa_bam_bai_ch
}

workflow filter_bam {

    take : 
        data

    main : 

        FILTER_BAM( data )
        processed_bam_ch = FILTER_BAM.out.filtered_bams_ch
        bam_to_bw_input_ch = FILTER_BAM.out.bam_to_bw_input

    emit : 
        filter_bam_bai = bam_to_bw_input_ch

}

workflow deepvariant {

    take : 
        data 

    main :
        DEEPVARIANT_CALL_VARIANTS( data, params.genome )

    emit : 
        vcfs = DEEPVARIANT_CALL_VARIANTS.out.deepvariant_vcf_ch

}
