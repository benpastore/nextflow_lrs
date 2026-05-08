#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { LONGQC_RAW } from '../modules/longqc/main.nf'
include { NANOPLOT_RAW } from '../modules/nanoplot/main.nf'
include { CHOPPER } from '../modules/chopper/main.nf'
include { MINIMAP2_INDEX } from '../modules/minimap2/main.nf'
include { MINIMAP2_ALIGN } from '../modules/minimap2/main.nf'
include { ALIGN_HAPLOTYPES_TO_REFERENCE } from '../minimap2/main.nf'
include { BAM_TO_BW } from '../modules/deeptools/main.nf'
include { INDEX_REFERENCE } from '../modules/samtools/main.nf'
include { FLYE_ASSEMBLE } from '../modules/flye/main.nf'
include { INDEX_FLYE_ASSEMBLY } from '../modules/samtools/main.nf'
include { HAPDUP } from '../modules/hapdup/main.nf'
include { CLAIR3 } from '../modules/clair3/main.nf'
include { WHATSHAP_PHASE } from '../modules/whatshap/main.nf'
include { WHATSHAP_HAPLOTAG } from '../modules/whatshap/main.nf'
include { SNIFFLES } from '../modules/sniffles/main.nf'
include { DIPCALL } from '../modules/dipcall/main.nf'
include { MAP_READS_TO_ASSEMBLY } from '../modules/minimap2/main.nf'

workflow index_reference { 
    take : 
        data 
    
    main : 
        INDEX_REFERENCE( data )
    
    emit : 
        reference_indexed_ch = INDEX_REFERENCE.out.reference_indexed_ch
}

workflow longqc_raw {
    take : 
        data
    
    main : 
        LONGQC_RAW(data)
}

workflow nanoplot_raw {
    take : 
        data
    
    main : 
        NANOPLOT_RAW(data)
}

// Reusable QC bundle: runs LongQC + NanoPlot on any reads channel
workflow qc_bundle_raw {
    take:
        data

    main:
        //LONGQC_RAW(data)
        NANOPLOT_RAW(data)
}

workflow qc_bundle_post {
    take:
        data

    main:
        //LONGQC_RAW(data)
        NANOPLOT_RAW(data)
}

workflow chopper { 
    take : 
        data
    
    main : 
        CHOPPER(data)
    
    emit : 
        reads = CHOPPER.out.fq
}

workflow minimap2 {

    take : 
        data

    main :
        genome_fasta = file("${params.genome}")
        genome_name = "${genome_fasta.baseName}"
        index_dir = "${params.index}/minimap2"
        index_path = "${params.index}/minimap2/${genome_name}.mmi"

        index_exists = file(index_path).exists()

        if ( index_exists == false || params.force_build_mmi){
            build = true
        } else {
            build = false
        }

        if (build == true){
            MINIMAP2_INDEX( params.genome, index_dir, index_path )
            index = MINIMAP2_INDEX.out.index_ch
        } else {
            index = index_path 
        }

        MINIMAP2_ALIGN( index, data )
    
    emit :
        bams = MINIMAP2_ALIGN.out.bam_ch
}

workflow bam_coverage { 
    take : 
        data
    
    main : 
        BAM_TO_BW( data )
    
    emit : 
        bws = BAM_TO_BW.out.bws_ch
}

workflow medaka { 
    take : 
        data
        reference
    
    main : 
        MEDAKA( data, reference )
}

workflow flye { 
    take : 
        data
        genome
    
    main : 
        FLYE_ASSEMBLE( data, genome )
    
    emit : 
        flye_assembly = FLYE_ASSEMBLE.out.flye_assembly
}

workflow polish_assembly_with_illumina (
    #################
)

workflow map_reads_to_flye_assembly {
    take : 
        data 
    
    main : 
        MAP_READS_TO_ASSEMBLY( data )
    
    emit : 
        reads_mapped_to_assembly = MAP_READS_TO_ASSEMBLY.out.aligned_reads_to_assembly
        //tuple val(sampleID), path(asm_fa), path(asm_fai), path("*.bam"), path("*.bai")
}

workflow hapdup_phase {
    take : 
        data
    
    main :
        HAPDUP( data )
    
    emit : 
        hapdup_ch = HAPDUP.out.hapdup_ch
}

workflow clair3 { 
    take : 
        data 
    main : 
        CLAIR3( params.genome, data )
    emit : 
        clair3_ch = CLAIR3.out.clair3_ch
}

workflow deepvariant {

    take : data 

    main :
        DEEPVARIANT_CALL_VARIANTS( data, params.genome )

    emit : 
        deepvariant_ch = DEEPVARIANT_CALL_VARIANTS.out.vcf

}

worfklow whatshap_phase { 
    take : 
        ref 
        data
    main : 
        WHATSHAP_PHASE( ref, data ) 
    emit : 
        whatshap_phase_ch = WHATSHAP_PHASE.out.whatshap_phase_ch
}

workflow whatshap_haplotag { 

    take : 
        ref
        data
    main : 
        WHATSHAP_HAPLOTAG( ref, data )
    emit : 
        whatshap_haplotag = WHATSHAP_HAPLOTAG.out.whathap_haplotag_ch

}

workflow sniffles { 

    take : 
        data 
    main : 
        SNIFFLES( data )
    emit : 
        sniffles = SNIFFLES.out.sniffles_ch

}

workflow align_haplotype_to_reference { 

    take : 
        ref
        data 
    
    main : 
        ALIGN_HAPLOTYPES_TO_REFERENCE( ref, data )
    
    emit : 
        aligned_haplotypes = ALIGN_HAPLOTYPES_TO_REFERENCE.out.aligned_haplotypes
    
}

workflow dipcall { 

    take : 
        ref
        data 
    
    main : 
        DIPCALL( ref, data)


}