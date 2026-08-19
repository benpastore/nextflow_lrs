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
    ========================================================================================
                            haplotype_aware long-read sequencing pipeline
    ========================================================================================
    Usage:
    The typical command for running the pipeline is as follows:

        nextflow main.nf -profile cluster --design design.csv --genome ref.fa --results results/

    or, starting from raw pod5 instead of an already-basecalled bam/fastq:

        nextflow main.nf -profile cluster --pod5_design pod5_design.csv --genome ref.fa --results results/

    Mandatory arguments:
    --results [dir]                     Path to results directory. All per-tool outputs are published under this.
    --genome [file]                     Reference genome fasta (T2T-CHM13 by default -- see nextflow.config).
    -profile [str]                      Configuration profile to use. Can use local / cluster.

    One of the following two is required to specify sample input:
    --design [file]                     Tab-separated file with columns SAMPLE, ONT, ONT_UNAL_BAM -- one row per
                                         sample (already-basecalled fastq + unaligned bam). Multiple rows per
                                         sample are merged (multi-lane/multi-run support).
    --pod5_design [file]                Alternative entry point: comma-separated csv (header: sample,pod5) mapping
                                         each sampleID to a directory of pod5 files. The pipeline basecalls with
                                         dorado (see --dorado_basecall_model) instead of reading --design.
                                         Mutually exclusive with --design -- if both are set, --pod5_design wins.

    Basecalling / trimming:
    --dorado_basecall_model [str]       dorado basecaller model string, used only with --pod5_design.
                                         (default: 'sup,5mCG_5hmCG' -- the 5mCG_5hmCG tag enables methylation
                                         calling; see the methylation section below)
    --default_dorado_seq_kit [str]      Sequencing kit passed to dorado trim/basecaller (default: 'SQK-LSK114')
    --dorado_trim [bool]                Run dorado trim before downstream processing (default: true)
    --dorado [bool]                     Run dorado polish during assembly (default: true)

    Assembly / alignment:
    --assemble_genome [bool]            Run the hifiasm assembly branch (dipcall/hapdiff/etc.) (default: true)
    --alignment_based_variant_calling [bool]
                                         Run the minimap2/clair3/longphase alignment branch (default: true)

    Error correction (HERRO -- GPU deep-learning ONT self-correction, https://github.com/lbcb-sci/herro):
    --herro_model [file]                Path to HERRO model weights, resolved inside the herro_inference
                                         container. Defaults to the R10.4.1 model baked into that image
                                         at build time (see docker/herro/Dockerfile.inference) -- no
                                         separate download needed. Runs on every sample, no per-sample
                                         gating (herro has no notion of matched Illumina data).
    --herro_gpu_ids [str]                GPU device IDs passed to `herro inference -d` (default: '0')
    --herro_batch_size [int]            `herro inference -b` (default: 64 -- 64 for 40GB VRAM, 128 for 80GB)
    --herro_preprocess_parts [int]      Chunk count for herro's preprocess.sh, bounds memory on large read
                                         sets (default: 4)

    Trio / haplotype family structure (hifiasm trio-binned assembly + WhatsHap pedigree phasing):
    --family_json [file]                JSON mapping each child sampleID to its parents:
                                         {"child": {"father": "...", "mother": "..."}}. Parents need no
                                         Illumina data -- they're identified purely by this file and must just be
                                         ordinary rows in the same --design/--pod5_design as children. (default:
                                         false, disabled)
    --run_trio_analysis [bool]          Master on/off switch for all trio-specific behavior. Requires
                                         --family_json to also be set. Leave false for singular-genome runs.
                                         (default: false)

    Paraphase (SMA/SMN1-SMN2 and other segmental-duplication region reconstruction):
    --paraphase_genome_build [str]      Reference build passed to paraphase's --genome (default: 'chm13')
    --paraphase_genes [str]             Comma-separated region names to restrict paraphase to (default: '',
                                         all regions supported for the chosen build)
    --paraphase_args [str]              Extra raw args appended to the paraphase command line

    Methylation:
    Enabled automatically whenever the input bam carries MM/ML tags (dorado's default mod-calling model, or
    --dorado_basecall_model with a mod-calling tag for --pod5_design runs). modkit pileup runs on the final
    haplotagged bam, split by haplotype (HP tag), producing per-haplotype bedMethyl under
    <results>/methylation/modkit. No separate flag needed to turn it on.

    AlphaGenome variant-effect annotation (scaffold only -- see bin/run_alphagenome.py TODOs):
    --run_alphagenome [bool]            (default: false)
    --alphagenome_weights [dir]         Required if --run_alphagenome is true
    --alphagenome_args [str]            Extra raw args appended to the alphagenome command line

    Misc:
    --publish_mode [str]                Nextflow publishDir mode used by most processes (default: 'copy')
    --outprefix [str]                   Output prefix (default: 'ONT_ANALYSIS' if unset)
    --debug [bool]                      Verbose logging / channel .view() calls (default: true)
    --help                              Show this message
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
include { hifasim_trio } from '../../subworkflows/ont.nf'
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
include { paraphase } from '../../subworkflows/ont.nf'
include { modkit } from '../../subworkflows/ont.nf'
include { dipcall } from '../../subworkflows/ont.nf'
include { hapdiff } from '../../subworkflows/ont.nf'
include { genome_stats } from '../../subworkflows/ont.nf'
include { chrom_coverage } from '../../subworkflows/ont.nf'
include { consolidate_variants } from '../../subworkflows/ont.nf'
include { alphagenome } from '../../subworkflows/ont.nf'

include { herro_correction } from '../../subworkflows/herro.nf'
include { trio_family } from '../../subworkflows/trio.nf'
include { WHATSHAP_PHASE_TRIO } from '../../modules/whatshap/main.nf'
include { pod5_basecall } from '../../subworkflows/pod5.nf'

// to add
include { DORADO_TRIM } from '../../modules/dorado/main.nf'
include { SAMTOOLS_CONVERT_BAM_TO_FASTQ } from '../../modules/samtools/main.nf'
include { SAMTOOLS_FASTQ_TO_BAM } from '../../modules/samtools/main.nf'
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
    if (params.pod5_design) {
        ch_pod5_design = file(params.pod5_design, checkIfExists: true)
    } else if (params.design) {
        ch_design = file(params.design, checkIfExists: true)
    } else {
        exit 1, 'Neither --design nor --pod5_design specified!'
    }
    if (params.genome)    { ch_genome = file(params.genome, checkIfExists: true) } else { exit 1, 'Genome fasta not specified!' }
    if (params.outprefix) { ; } else {'Outprefix not specified! Defaulting to ONT_ANALYSIS'; params.outprefix = 'ONT_ANALYSIS' }

    // upstream of everything make sure minimap2 index is built for reference
    minimap2_index( params.genome, params.index ) // THIS GENERATES THE MINIMAP2 ALIGNMENT INDEX
    reference_index = minimap2_index.out.index // THIS GENERATES THE MINIMAP2 ALIGNMENT INDEX

    INDEX_REFERENCE( params.genome ) // THIS IS THE SAMTOOLS INDEX .fai GENERATION NOT FOR ALIGNMENT
    samtools_fai_index = INDEX_REFERENCE.out.ref_indexed_ch // THIS IS THE SAMTOOLS INDEX .fai GENERATION NOT FOR ALIGNMENT

    ////////////////// parse design (or basecall from pod5, if pod5_design is set)
    if (params.pod5_design) {
        pod5_basecall( params.pod5_design )
        ont_reads = pod5_basecall.out.ont
    } else {
        parse_design( params.design )
        ont_reads = parse_design.out.ont
    }

    if (params.debug) { 
        // for debugging view ont reads
        ont_reads.view() 
    }

    // trimming upstream of assembly and alignment
    println("Dorado trim sequencing kit set to ${params.default_dorado_seq_kit}")
    DORADO_TRIM( ont_reads )

    println("Samtools convert to fastq")
    SAMTOOLS_CONVERT_BAM_TO_FASTQ( DORADO_TRIM.out.dorado_trim_output_ch )

    //ont_unal_bam = DORADO_TRIM.out.dorado_trim_output_ch
    ont_ch = SAMTOOLS_CONVERT_BAM_TO_FASTQ.out.samtool_convert_unalbam2fastq

    /////////////////// prelim QC on nanopore reads (always on raw reads, before HERRO correction)
    NANOPLOT_RAW( ont_ch )

    ////////////////// GPU DEEP-LEARNING ERROR CORRECTION (HERRO) /////////////////////
    // every sample gets herro-corrected (no per-sample gating -- herro is ONT-only
    // self-correction, unlike the old Illumina-hybrid HERO tool). Feeds assembly and
    // alignment below.
    herro_correction( ont_ch )
    ont_ch_corrected = herro_correction.out.ont_ch

    ////////////////// TRIO / HAPLOTYPE FAMILY STRUCTURE /////////////////////
    // params.run_trio_analysis is the single on/off switch for all
    // trio-specific behavior (hifiasm trio-binning + WhatsHap pedigree
    // phasing below). Flip it off for singular-genome runs without having
    // to also unset/remove family_json. When off (or family_json unset),
    // trio_family is never called and every child just takes the regular
    // non-trio assembly path below.
    //
    // Parents have no Illumina data -- they're just ordinary rows in the
    // same ONT design.csv/pod5_design as children, identified purely via
    // family_json. trio_family splits them out into non_parent_ont_ch so
    // they never enter assembly/alignment/phasing/paraphase/modkit below as
    // full samples; only ont_ch_corrected gets swapped, everything
    // downstream of it is unchanged.
    run_trio = params.run_trio_analysis && params.family_json
    if (run_trio) {
        // computed once, reused at both the assembly gate below and the
        // phasing gate further down -- avoids re-parsing family_json or
        // re-running yak/clair3 for the same parents twice
        trio_family( params.family_json, ont_ch_corrected, reference_index, samtools_fai_index )
        ont_ch_corrected = trio_family.out.non_parent_ont_ch
    }

    ////////////////// ASSEMBLY SECTION /////////////////////
    // default (no-assembly) placeholders so the variant consolidation join below
    // always has a channel to join against, even when assembly is skipped
    dipcall_vcf_ch = Channel.empty()
    hapdiff_vcf_ch = Channel.empty()

    if (params.assemble_genome) {
        //Hifasim: assemble ont reads
        if (params.debug) { println("ASSEMBLE GENOME WITH HIFASIM") }

        if (run_trio) {
            // children with a fully-resolved trio (both parents' yak dbs
            // available) get hifiasm trio-binning; everyone else gets the
            // regular non-trio assembly. Merged back into one channel shaped
            // like hifasim.out.hifasim_asm so nothing further downstream
            // (gfatools/dorado/dipcall/hapdiff) needs to know which path a
            // given sample took.
            assembly_gate = ont_ch_corrected
                .join( trio_family.out.trio_assembly, remainder: true )
                .branch {
                    has_trio: it[3] != null
                    no_trio: true
                }

            hifasim( assembly_gate.no_trio.map { sid, bam, fastq, pat_yak, mat_yak -> tuple(sid, bam, fastq) } )

            hifasim_trio( assembly_gate.has_trio.map { sid, bam, fastq, pat_yak, mat_yak -> tuple(sid, bam, fastq, pat_yak, mat_yak) } )

            hifasim_asm_ch = hifasim.out.hifasim_asm.mix( hifasim_trio.out.hifasim_trio_asm )
        } else {
            // singular-genome mode: no trio structure, everyone takes the
            // regular non-trio assembly path
            hifasim( ont_ch_corrected )
            hifasim_asm_ch = hifasim.out.hifasim_asm
        }

        // use jasmine to combine variant calls from multiple callers
        // use busco to assess quality of assembly

        //convert hifasim to fasta
        gfatools( hifasim_asm_ch )

        // include dorado polish
        if (params.dorado) {

            dorado( gfatools.out.fasta_asm )
            genome_asm_ch = dorado.out.dorado_output_ch

        } else {
            genome_asm_ch = gfatools.out.haplotype_asm
                .map { sampleID, bam, reads, hap1, hap2 -> tuple(sampleID, reads, hap1, hap2) }
        }

        // aggregate read/base counts and assembled genome size for QC
        genome_stats( genome_asm_ch )
        
        // map the assembly to the reference
        minimap2_map_asm_to_ref( reference_index, genome_asm_ch)

        // per-chromosome average coverage of the assembly-vs-reference alignment
        chrom_coverage( minimap2_map_asm_to_ref.out.bams )

        //run dip call & hapdiff
        // index the reference genome quickly for dipcall 
        dipcall( genome_asm_ch,  samtools_fai_index )
        hapdiff( genome_asm_ch,  samtools_fai_index )

        // assembly-based variant calls, for the master variant table below
        dipcall_vcf_ch = dipcall.out.dipcall
            .map { sampleID, vcf, tbi -> tuple(sampleID, vcf) }
        hapdiff_vcf_ch = hapdiff.out.hapdiff_ch

    }
    ////////////////// VARIANT CALLING SECTION -- From alignment & haplotype phasing /////////////////////
    params.alignment_based_variant_calling = true
    if (params.alignment_based_variant_calling) { 

        // MAP READS TO REFERENCE GENOME
        //minimap2( reference_index, ont_reads )
        minimap2(reference_index, ont_ch_corrected.map { sid, unal_bam, fastq -> tuple(sid, fastq) })
        
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

        // paraphase (SMA/SMN1-SMN2 region reconstruction from segmental duplications)
        paraphase( haplo_ch_v2, samtools_fai_index )

        // methylation pileup (5mCG/5hmCG), split by haplotype via longphase_sv's HP tags
        modkit( haplo_ch_v2, samtools_fai_index )

        ////////////////// CONSOLIDATE VARIANTS SECTION /////////////////////
        // master per-sample VCF: alignment-based (clair3 small variants,
        // phased+SV-aware via longphase_sv) + sniffles (SV) callers, joined
        // with assembly-based (dipcall, hapdiff) callers when available.
        // Each caller keeps its own (phase-aware) genotype column so
        // haplotype info is preserved; INFO/CALLERS + INFO/NUM_CALLERS
        // record which/how many callers detected each variant.
        clair3_longphase_vcf_ch = longphase_sv.out.longphase_sv_vcf_ch
            .map { sampleID, vcf, tbi -> tuple(sampleID, vcf) }

        sniffles_vcf_ch = sniffles.out.sniffles_ch
            .map { sampleID, vcf, tbi -> tuple(sampleID, vcf) }

        consolidate_input_ch = clair3_longphase_vcf_ch
            .join(sniffles_vcf_ch)
            .join(dipcall_vcf_ch, remainder: true)
            .join(hapdiff_vcf_ch, remainder: true)
            .map { sampleID, longphase_vcf, sniffles_vcf, dipcall_vcf, hapdiff_vcf ->
                tuple(
                    sampleID,
                    longphase_vcf,
                    sniffles_vcf,
                    dipcall_vcf ?: file("NO_FILE"),
                    hapdiff_vcf ?: file("NO_FILE")
                )
            }

        consolidate_variants( consolidate_input_ch )

        ////////////////// TRIO-AWARE PEDIGREE PHASING (children only) /////////////////////
        // Runs after variant consolidation (children-only step; parents never
        // appear in longphase_sv/haplo_ch_v2 at all, since they're never run
        // as full samples -- see trio_family). Additive refinement of the
        // child's longphase_sv-phased small variants using Mendelian
        // constraints from both parents' genotype VCFs (WhatsHap --ped).
        // Phases the same clair3/longphase small-variant VCF consolidate_variants
        // consumed above, not the master.vcf itself -- that file's sample
        // columns are caller names, not sample IDs, and it mixes in symbolic
        // SV records neither of which WhatsHap can phase. This step doesn't
        // feed back into consolidate_variants/longphase_sv/straglr/spectre --
        // it's purely an extra, separately-published output.
        if (run_trio) {
            trio_phasing_input = longphase_sv.out.longphase_sv_vcf_ch
                .join( haplo_ch_v2 )
                .join( trio_family.out.trio_phasing )
                .map { sampleID, vcf, tbi, bam, bai, father, mother, father_vcf, mother_vcf ->
                    tuple(sampleID, father, mother, vcf, bam, bai, father_vcf, mother_vcf)
                }

            WHATSHAP_PHASE_TRIO( samtools_fai_index, trio_phasing_input )
        }

        ////////////////// ALPHAGENOME VARIANT-EFFECT ANNOTATION /////////////////////
        // scaffold only -- disabled by default until bin/run_alphagenome.py's
        // TODOs are filled in against your local weights directory
        if (params.run_alphagenome) {
            if (!params.alphagenome_weights) { exit 1, 'params.alphagenome_weights not set!' }

            alphagenome_input_ch = consolidate_variants.out.consolidated_ch

            alphagenome( alphagenome_input_ch, samtools_fai_index )
        }

    }

}