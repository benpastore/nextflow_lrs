#!/usr/bin/env nextflow

include { ATRIA } from '../modules/atria/main.nf'
include { BWA_INDEX } from '../modules/bwa/main.nf'
include { BWA_MEM } from '../modules/bwa/main.nf'
include { FILTER_BAM } from '../modules/samtools/main.nf'
include { INDEX_BAM } from '../modules/samtools/main.nf'
include { BAM_TO_BW } from '../modules/deeptools/main.nf'
include { DEEPVARIANT_CALL_VARIANTS } from '../modules/deepvariant/main.nf'

nextflow.enable.dsl=2

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
