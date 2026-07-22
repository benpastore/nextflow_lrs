#!/usr/bin/env nextflow

/*
========================================================================================
                         long read sequencing pipeline
========================================================================================
Ben Pastore
pastore.28@osu.edu
----------------------------------------------------------------------------------------
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
params.index = "${params.base}/index"
/*
////////////////////////////////////////////////////////////////////
Enable dls2 language --> import modules
////////////////////////////////////////////////////////////////////
*/
nextflow.enable.dsl=2

/* include modules */
include { DESIGN_INPUT } from '../../modules/parse_design/main.nf'
include { NANOPLOT_RAW } from '../../modules/nanoplot/main.nf'
include { INDEX_REFERENCE } from '../../modules/samtools/main.nf'

/* include workflows */
///////////// ILLUMINA ////////////
include { atria } from '../../subworkflows/illumina.nf'
include { bwa_align } from '../../subworkflows/illumina.nf'
include { filter_bam } from '../../subworkflows/illumina.nf'
include { deepvariant } from '../../subworkflows/illumina.nf'

/////////////// ONT //////////////////
include { hifasim } from '../../subworkflows/ont.nf'
include { gfatools } from '../../subworkflows/ont.nf'
include { dorado } from '../../subworkflows/ont.nf'

include { minimap2 } from '../../subworkflows/ont.nf'
include { minimap2_index } from '../../subworkflows/ont.nf'
include { minimap2_map_asm_to_ref } from '../../subworkflows/ont.nf'

include { longphase } from '../../subworkflows/ont.nf'
include { longphase_sv } from '../../subworkflows/ont.nf'

include { clair3 } from '../../subworkflows/ont.nf'
include { sniffles } from '../../subworkflows/ont.nf'
include { spectre } from '../../subworkflows/ont.nf'
include { straglr } from '../../subworkflows/ont.nf'
include { dipcall } from '../../subworkflows/ont.nf'
include { hapdiff } from '../../subworkflows/ont.nf'

// to add 
include { DORADO_TRIM } from '../../modules/dorado/main.nf'
include { SAMTOOLS_CONVERT_BAM_TO_FASTQ } from '../../modules/samtools/main.nf'
include { PORECHOP } from '../../modules/porechop/main.nf'
include { MERGE_INPUTS } from '../../modules/samtools/main.nf'

workflow parse_design {
    take : 
        data

    main : 
        DESIGN_INPUT( data )

        ont = DESIGN_INPUT
            .out
            .fastq_ch
            .splitCsv( header: ['sample', 'ont', 'unal_bam'], sep: ",", skip: 1)
            .map{ row -> [ row.sample, file(row.ont), file(row.unal_bam) ] }
            .groupTuple()
        
        ont.view()

        MERGE_INPUTS( ont )
        
        ont_ch = MERGE_INPUTS.out.merged_inputs_ch

        /*
        ont = DESIGN_INPUT
            .out
            .fastq_ch
            .splitCsv( header: ['sample', 'ont', 'unal_bam', 'r1', 'r2'], sep: ",", skip: 1)
            .map{ row -> [ row.sample, row.ont ] }
        
        ont_unal_bam = DESIGN_INPUT
            .out
            .fastq_ch
            .splitCsv( header: ['sample', 'ont', 'unal_bam', 'r1', 'r2'], sep: ",", skip: 1)
            .map{ row -> [ row.sample, row.unal_bam ] }

        ont_illumina = DESIGN_INPUT
            .out
            .fastq_ch
            .splitCsv( header: ['sample', 'ont', 'unal_bam', 'r1', 'r2'], sep: ",", skip: 1)
            .map{ row -> [ row.sample, row.ont, row.r1, row.r2 ] }

        illumina = DESIGN_INPUT
            .out
            .fastq_ch
            .splitCsv( header: ['sample', 'ont', 'unal_bam', 'r1', 'r2'], sep: ",", skip: 1)
            .map{ row -> [ row.sample, row.r1, row.r2 ] }

        //replicates_ch = DESIGN_INPUT
        //    .out
        //    .condition_ch
        //    .splitCsv(header:false, skip:1)
        */
        
    emit : 
        //ont = ont
        //ont_illumina = ont_illumina
        //illumina = illumina
        ont = ont_ch
        //replicates = replicates_ch
}

/*
//////////////////////// MAIN WORKFLOW ////////////////////////
*/

workflow { 

    /*
    To Do: 
    - Dorado2 base calling with methylation detection
    */

    /*
    ////////////////////////////////////////////////////////////////////
    Validate mandatory inputs (design, genome, junctions, results, outprefix)
    ////////////////////////////////////////////////////////////////////
    */
    if (params.design)    { ch_design = file(params.design, checkIfExists: true) } else { exit 1, 'Design file not specified!' }
    if (params.genome)    { ch_genome = file(params.genome, checkIfExists: true) } else { exit 1, 'Genome fasta not specified!' }
    if (params.outprefix) { ; } else {'Outprefix not specified! Defaulting to ONT_ANALYSIS'; params.outprefix = 'ONT_ANALYSIS' }

    // upstream of everything make sure minimap2 index is built for reference 
    minimap2_index( params.genome, params.index ) // THIS GENERATES THE MINIMAP2 ALIGNMENT INDEX
    reference_index = minimap2_index.out.index // THIS GENERATES THE MINIMAP2 ALIGNMENT INDEX

    INDEX_REFERENCE( params.genome ) // THIS IS THE SAMTOOLS INDEX .fai GENERATION NOT FOR ALIGNMENT
    samtools_fai_index = INDEX_REFERENCE.out.ref_indexed_ch // THIS IS THE SAMTOOLS INDEX .fai GENERATION NOT FOR ALIGNMENT

    ////////////////// parse design
    parse_design( params.design )
    ont_reads = parse_design.out.ont

    //ont_unal_bam = parse_design.out.ont_unal_bam
    //ont_illumina = parse_design.out.ont_illumina
    //illumina_reads = parse_design.out.illumina
    
    if (params.debug) { 
        // for debugging view ont reads
        ont_reads.view() 
    }

    // dorado2 basecall
    // input: sample, /path/to/pod5
    // output: bam / fastq
    // combine sample bam/fastq

    // trimming upstream of assembly and alignment
    println("Dorado trim sequencing kit set to ${params.default_dorado_seq_kit}")
    DORADO_TRIM( ont_reads )

    println("Samtools convert to fastq")
    SAMTOOLS_CONVERT_BAM_TO_FASTQ( DORADO_TRIM.out.dorado_trim_output_ch )

    //println("Porechop")
    PORECHOP( SAMTOOLS_CONVERT_BAM_TO_FASTQ.out.samtool_convert_unalbam2fastq )

    //ont_unal_bam = DORADO_TRIM.out.dorado_trim_output_ch
    ont_ch = SAMTOOLS_CONVERT_BAM_TO_FASTQ.out.samtool_convert_unalbam2fastq
    ont_porechop = PORECHOP.out.porechop_ch

    /*
    /////////////////// ILLUMINA SECTION **OPTIONAL** /////////////
    if (params.use_illumina) {

        if (params.debug) { println("RUNNING ILLUMINA PIPELINE") }

        // atria
        atria( illumina_reads )

        // bwa 
        bwa_align( atria.out.fqs )

        // filter alignments
        filter_bam( bwa_align.out.bams )

        // deepvariant
        deepvariant( filter_bam.out.filter_bam_bai )

    } else { 
        if (params.debug) { println("SKIPPING ILLUMINA") }
    }
    */
        
    /////////////////// prelim QC on nanopore reads 
    NANOPLOT_RAW( ont_ch )

    ////////////////// ASSEMBLY SECTION /////////////////////
    if (params.assemble_genome) { 
        //Hifasim: assemble ont reads  
        if (params.debug) { println("ASSEMBLE GENOME WITH HIFASIM") }
        hifasim( ont_ch )

        // Include an if else here to use verkko instead (may not need to use gfatools, depending on verkko output)
        // assembly QC & polishing with medaka
        // # contigs, % completeness, contig size 
        // go directly into dipcall from hifiasm 
        // minimira --> remove ligation based artifacts
        // use jasmine to combine variant calls from multiple callers
        // use busco to assess quality of assembly

        //convert hifasim to fasta 
        gfatools( hifasim.out.hifasim_asm )

        // include dorado polish
        if (params.dorado) {

            // dorado needs unaligned bam to work, so merge unaligned bam into gfatools.out
            //dorado_input_ch = gfatools
            //    .out
            //    .fasta_asm
            //    .join(ont_unal_bam)
            //    .map { sampleID, reads, hapfasta, hapfai, unal_bam  ->
            //        tuple( sampleID, reads, hapfasta, hapfai, unal_bam )
            //     }
            //dorado_input_ch.view { x -> "DORADO INPUT: $x" }

            dorado( gfatools.out.fasta_asm )
            genome_asm_ch = dorado.out.dorado_output_ch

        } else {
            genome_asm_ch = gfatools.out.haplotype_asm
                .map { sampleID, bam, reads, hap1, hap2 -> tuple(sampleID, reads, hap1, hap2) }
        }

        // map the asmbley to the reference
        minimap2_map_asm_to_ref( reference_index, genome_asm_ch)

        //run dip call & hapdiff
        // index the reference genome quickly for dipcall 
        dipcall( genome_asm_ch,  samtools_fai_index )
        hapdiff( genome_asm_ch,  samtools_fai_index )

    }
    ////////////////// VARIANT CALLING SECTION -- From alignment & haplotype phasing /////////////////////
    params.alignment_based_variant_calling = true
    if (params.alignment_based_variant_calling) { 

        // MAP READS TO REFERENCE GENOME
        //minimap2( reference_index, ont_reads )
        minimap2(reference_index, ont_porechop)
        
        // clair3 call variants 
        clair3( samtools_fai_index, minimap2.out.bams )

        variants_ch = clair3.out.clair3_ch

        vcf_bam_ch = variants_ch
            .join(minimap2.out.bams)
            .map { sampleID, vcf, tbi, bam, bai ->
                tuple(sampleID, vcf, tbi, bam, bai)
            }

        // phase using variants from clair3
        longphase( samtools_fai_index, vcf_bam_ch )
        haplo_ch = longphase.out.longphase_ch

        // sniffles (SV calling)
        sniffles( haplo_ch )

        // re-phase with long phase including sv
        // combine sniffles + clair3 -> do another long phase and pass that to stragglr and spectre
        clair3_sniffles_bam_ch = clair3.out.clair3_ch
            .join(sniffles.out.sniffles_ch)
            .join(minimap2.out.bams)
            .map { sampleID, clair3_vcf, clair3_tbi, sniffles_vcf, sniffles_tbi, bam, bai ->
                tuple(sampleID, clair3_vcf, clair3_tbi, sniffles_vcf, sniffles_tbi, bam, bai)
            }
        
        // longphase round2
        longphase_sv( samtools_fai_index, clair3_sniffles_bam_ch )
        haplo_ch_v2 = longphase_sv.out.longphase_sv_ch

        // straggler (expansion repeats)
        straglr( haplo_ch_v2, samtools_fai_index )

        // spectre CNV
        spectre( haplo_ch_v2, samtools_fai_index )

    }

}