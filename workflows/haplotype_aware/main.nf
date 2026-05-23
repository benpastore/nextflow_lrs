#!/usr/bin/env nextflow


/*
Discussion points for meeting

- hifiasm vs verkko 
- pepper, margin?
- quality assesment of assembly
/*

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
include { map_reads_to_assembly } from '../../subworkflows/ont.nf'
include { hapdup_phase } from '../../subworkflows/ont.nf'
include { INDEX_REFERENCE } from '../../modules/samtools/main.nf'
include { minimap2 } from '../../subworkflows/ont.nf'
include { clair3 } from '../../subworkflows/ont.nf'
include { whatshap_phase } from '../../subworkflows/ont.nf'
include { whatshap_haplotag } from '../../subworkflows/ont.nf'
include { sniffles } from '../../subworkflows/ont.nf'
include { spectre } from '../../subworkflows/ont.nf'
include { straglr } from '../../subworkflows/ont.nf'
include { dipcall } from '../../subworkflows/ont.nf'
include { NANOPLOT_RAW } from '../../modules/nanoplot/main.nf'

workflow parse_design {
    take : 
        data

    main : 
        DESIGN_INPUT( data )

        ont = DESIGN_INPUT
            .out
            .fastq_ch
            .splitCsv( header: ['sample', 'ont', 'r1', 'r2'], sep: ",", skip: 1)
            .map{ row -> [ row.sample, row.ont ] }

        ont_illumina = DESIGN_INPUT
            .out
            .fastq_ch
            .splitCsv( header: ['sample', 'ont', 'r1', 'r2'], sep: ",", skip: 1)
            .map{ row -> [ row.sample, row.ont, row.r1, row.r2 ] }

        illumina = DESIGN_INPUT
            .out
            .fastq_ch
            .splitCsv( header: ['sample', 'ont', 'r1', 'r2'], sep: ",", skip: 1)
            .map{ row -> [ row.sample, row.r1, row.r2 ] }

        //replicates_ch = DESIGN_INPUT
        //    .out
        //    .condition_ch
        //    .splitCsv(header:false, skip:1)
        
    emit : 
        ont = ont
        ont_illumina = ont_illumina
        illumina = illumina
        //replicates = replicates_ch
}

/*
//////////////////////// MAIN WORKFLOW ////////////////////////
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

    ////////////////// parse design

        parse_design( params.design )
        ont_reads = parse_design.out.ont
        ont_illumina = parse_design.out.ont_illumina
        illumina_reads = parse_design.out.illumina
        
        if (params.debug) { 
            // for debugging view ont reads
            ont_reads.view() 
        }

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
        
    /////////////////// prelim QC on nanopore reads 
        NANOPLOT_RAW( ont_reads )

    ////////////////// ASSEMBLY SECTION /////////////////////
        //1. Hifasim: assemble ont reads  
        if (params.debug) { println("ASSEMBLE GENOME WITH HIFASIM") }
        hifasim( ont_reads ) // check hifiasm version
        hifasim.out.hifasim_asm.view()

        // Include an if else here to use verkko instead (may not need to use gfatools, depending on verkko output)
        // assembly QC & polishing with medaka
        // # contigs, % completeness, contig size 
        // go directly into dipcall from hifiasm 

        // minimira --> remove ligation based artifacts
        // use jasmine to combine variant calls from multiple callers
        // use busco to assess quality of assembly

        //2. convert hifasim to fasta 
        //if (params.debug) { println("CONVERT HIFASIM ASM TO FASTA") }
        gfatools( hifasim.out.hifasim_asm )

        // include dorado polish
        if (params.dorado) { 
            dorado( gfatools.out.fasta_asm )
            genome_asm_ch = dorado.out.dorado_output_ch
        } else { 
            genome_asm_ch = gfatools.out.haplotype_asm
        }

        //2. run dip call
        // index the reference genome quickly for dipcall 
        INDEX_REFERENCE( params.genome )
        dipcall( genome_asm_ch,  INDEX_REFERENCE.out.ref_indexed_ch )
        dipcall.out.dipcall.view()

    ////////////////// VARIANT CALLING SECTION -- From alignment & haplotype phasing /////////////////////
        // MAP READS TO REFERENCE GENOME
        minimap2( ont_reads )

        // clair3 call variants 
        clair3( minimap2.out.bams)

        // spectre copy number variant caller 
        spectre_input_ch = minimap2.out.bams
            .combine( INDEX_REFERENCE.out.ref_indexed_ch )
            .map { sampleID, bam, bai, ref_fa, ref_fai ->
                tuple(sampleID, bam, bai, ref_fa, ref_fai)
            }
        spectre( spectre_input_ch )

        // mix the variant outputs
        if ( params.use_illumina ) {
            clair3_vcf_ch = clair3.out.clair3_ch
                .map { sampleID, vcf, tbi ->
                    tuple(sampleID, 'clair3', vcf, tbi)
                }
            deepvariant_vcf_ch = deepvariant.out.vcfs
                .map { sampleID, vcf, tbi ->
                    tuple(sampleID, 'deepvariant', vcf, tbi)
                }

            variants_ch = clair3_vcf_ch.mix(deepvariant_vcf_ch)
        } else { 
            variants_ch = clair3.out.clair3_ch
        }

        vcf_bam_ch = variants_ch
            .join(minimap2.out.bams)
            .map { sampleID, caller, vcf, tbi, bam, bai ->
                tuple(sampleID, caller, vcf, tbi, bam, bai)
            }
   
        // whats haplotype phase on clair3, in future can also run dv -> merge with clair3 -> run haplotype phasing
        // can also run something called longphase ( use long phase )
        whatshap_phase( params.genome, vcf_bam_ch )

        // whatshap happlotag
        whatshap_haplotag( params.genome, whatshap_phase.out.whatshap_phase_ch )

        // sniffles (SV calling)
        sniffles( whatshap_haplotag.out.whatshap_haplotag )

        // straggler
        straglr( whatshap_haplotag.out.whatshap_haplotag )

    /*
    - https://github.com/epi2me-labs/wf-human-variation
     - longphase (could be used instead of whathap)
     - hapdiff (think this does variant calling from the de novo assembled genomes)
     */
}