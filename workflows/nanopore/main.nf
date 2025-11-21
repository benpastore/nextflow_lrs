#!/usr/bin/env nextflow

/*
========================================================================================
                         long read sequencing pipeline
========================================================================================
Ben Pastore
pastore.28@osu.edu
----------------------------------------------------------------------------------------


Softwares:
> longQC --> quality assesment 
> nanoplot --> quality assesment
> porechop --> remove adapters
> filtlong --> quality and length based filtering 
> nanofilt --> trimming ends and general filtering
> minimap2 --> alignment 
> DeepTools --> conversion of bam to bw file
> DeepVariant --> variant calling

*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:

    nextflow main.nf -profile cluster --design design.csv 

    Mandatory arguments:
    --design [file]                     Comma-separated file containing information about the samples in the experiment (see docs/usage.md) (Default: './design.csv')
    --results [file]                    Path to results directory
    -profile [str]                      Configuration profile to use. Can use local / cluster
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

params.bin = "${params.base}/../../bin"
params.index = "${params.base}/../../index"
/*
////////////////////////////////////////////////////////////////////
Enable dls2 language --> import modules
////////////////////////////////////////////////////////////////////
*/
nextflow.enable.dsl=2

include { DESIGN_INPUT } from '../../modules/parse_design/main.nf'
include { LONGQC_RAW } from '../../modules/longqc/main.nf'
include { NANOPLOT_RAW } from '../../modules/nanoplot/main.nf'
include { CHOPPER } from '../../modules/chopper/main.nf'
include { MINIMAP2_INDEX } from '../../modules/minimap2/main.nf'
include { MINIMAP2_ALIGN } from '../../modules/minimap2/main.nf'
include { BAM_TO_BW } from '../../modules/deeptools/main.nf'

workflow parse_design {

    take : 
        data

    main : 
        DESIGN_INPUT( data )

        reads_ch = DESIGN_INPUT
            .out
            .fastq_ch
            .splitCsv( header: ['sample', 'reads'], sep: ",", skip: 1)
            .map{ row -> [ row.sample, row.reads ] }

        replicates_ch = DESIGN_INPUT
            .out
            .condition_ch
            .splitCsv(header:false, skip:1)
        
    emit : 
        reads = reads_ch

        replicates = replicates_ch

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

/*
////////////////////////////////////////////////////////////f=////////
workflow
////////////////////////////////////////////////////////////////////
*/

workflow {

    /*
    ////////////////////////////////////////////////////////////////////
    Validate mandatory inputs (design, genome, junctions, results, outprefix)
    ////////////////////////////////////////////////////////////////////
    */
    if (params.design)    { ch_design = file(params.design, checkIfExists: true) } else { exit 1, 'Design file not specified!' }
    if (params.genome)    { ch_genome = file(params.genome, checkIfExists: true) } else { exit 1, 'Genome fasta not specified!' }
    if (params.outprefix) { ; } else {'Outprefix not specified! Defaulting to ONT_ANALYSIS'; params.outprefix = 'ONT_ANALYSIS' }

    // parse design 
    parse_design( params.design )

    // longqc 
    // nanoplot 
    parse_design.out.reads.view()
    qc_bundle_raw( parse_design.out.reads )

    // chopper 
    chopper( parse_design.out.reads )

    // longqc post 
    // nanoplot post
    qc_bundle_post( chopper.out.reads )

    // minimap2
    minimap2( chopper.out.reads )

    //bam to bw
    bam_coverage( minimap2.out.bams )
}




