#!/usr/bin/env nextflow

/*
========================================================================================
                         long read sequencing pipeline
========================================================================================
Ben Pastore
pastore.28@osu.edu

Flow (each stage is a numbered section in the workflow below):
  1. reads        parse design / basecall pod5 -> dorado trim -> QC -> HERRO correction
  2. family       split parents out (family.json); parents get alignment + Clair3 + Sniffles only
  3. alignment    minimap2 -> Clair3 -> longphase -> Sniffles -> longphase (SNV+SV)
  4. other calls  Straglr repeats, Spectre CNV, paraphase, methylation (optional)
  5. trio phasing WhatsHap pedigree phasing -> the child's phased small-variant VCF
                  (SVs: SV_TRIO_ORIGIN, since WhatsHap can't phase SVs)
  6. annotate     consolidate (WhatsHap small variants + Sniffles SVs), every variant
                  labelled with its parental haplotype (INFO/ORIGIN) -> snpEff -> population AF
  7. IBD          hap-ibd on the family-phased VCF: which parental copies each mild /
                  severe sib pair shares (IBD0/1/2)
  8. modifiers    per family, the mild child's genes with a homozygous variant, a
                  compound het or a non-genic (regulatory) variant that differs from
                  the severe sib -> mild / severe VCFs with IBD state + AF
  9. methylation  mild vs severe promoter methylation per parental copy (optional)
 10. recurrence   genes passing in >= modifier_min_mild mild children -> AlphaGenome /
                  AlphaMissense on every mild and severe variant in the gene
                  (two deduplicated VCFs, scored per side)

  With parents:
    WHATSHAP_PHASE_FAMILY ---> family VCF -> ORIENT_PHASE -> HAP_IBD -> IBD_SEGMENTS --+
         |                                                                            +-> MODIFIER_VARIANTS
         +-> child VCFs --+                                                           |    -> GENE_RECURRENCE
    SV_TRIO_ORIGIN -------+-> CONSOLIDATE (INFO/ORIGIN) -> SNPEFF -> ANNOTATE_AF -----+    -> AlphaGenome / AlphaMissense
    IBD_SEGMENTS -> METHYLATION_COMPARE
  (assembly: off by default, --assemble_genome true; see subworkflows/assembly.nf)
  Final mild vs severe outputs, <results>/10_modifier_analysis/:
    cnv/  repeat_expansion/  methylation/  variants/
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    ========================================================================================
                            haplotype_aware long-read sequencing pipeline
    ========================================================================================
    Usage:
        nextflow main.nf -profile cluster --design design.csv --genome ref.fa --results results/
        nextflow main.nf -profile cluster --pod5_design pod5_design.csv --genome ref.fa --results results/

    Mandatory arguments:
    --results [dir]                     Results directory
    --genome [file]                     Reference fasta (T2T-CHM13 by default -- see nextflow.config)
    -profile [str]                      local / cluster
    --design [file]                     Comma-separated: sample, ONT fastq, ONT unaligned bam (rows per sample are merged)
      or --pod5_design [file]           Comma-separated: sample,pod5 directory -- basecalls with dorado (wins if both set)

    Basecalling / trimming:
    --dorado_basecall_model [str]       dorado model, pod5 entry only (default 'sup,5mCG_5hmCG')
    --default_dorado_seq_kit [str]      Kit for dorado trim/basecaller (default 'SQK-LSK114')

    Error correction (HERRO):
    --herro_model [file]                HERRO R10.4.1 model weights
    --herro_gpu_ids [str] / --herro_batch_size [int] / --herro_preprocess_parts [int]

    Families:
    --family_json [file]                {"child": {"father": "...", "mother": "...", "phenotype": "mild|severe"}}
                                         Parents are ordinary rows in the same design file.
    --run_trio_analysis [bool]          Use the parents (family phasing, parent of origin, sibling
                                         haplotype sharing) (default: false)

    Sibling IBD with hap-ibd (<results>/07_haplotypes):
    --hapibd_jar [file]                 hap-ibd.jar (required with --run_trio_analysis)
    --hapibd_min_seed / --hapibd_min_extend / --hapibd_min_output [cM]   hap-ibd segment lengths (2.0 / 1.0 / 2.0)
    --hapibd_max_gap [bp]               hap-ibd max-gap (default 1000)
    --hapibd_min_gq [int]               drop child hets below this GQ before hap-ibd (default 20)
    --hapibd_bridge_cm [cM]             not_shared gaps shorter than this between shared segments are bridged (2.0)
    --hap_exclude_chroms [str]          Not scanned (default 'chrX,chrY,chrM')

    Modifier search (<results>/10_modifier_analysis/variants; needs --run_trio_analysis + --run_snpeff):
    --run_modifier_analysis [bool]      Run it (default: true)
    --modifier_sv_match_dist [bp]       SVs of the same type within this distance match (500)
    --modifier_min_mild [int]           Mild children a gene must pass in for AlphaGenome / AlphaMissense (2)
    --modifier_candidate_max_af [float] AlphaGenome / AlphaMissense input only: drop variants with every
                                         annotated AF >= this (default "" = no filter)
    --run_alphagenome / --run_alphagenome_api [bool]   AlphaGenome local / hosted API (default: false)
    --alphamissense_tsv [file]          AlphaMissense table; set = run AlphaMissense

    Optional branches:
    --assemble_genome [bool]            hifiasm assembly -> dipcall/hapdiff (default: false)
    --run_methylation [bool]            Raw-read alignment + modkit per haplotype (default: false); with
                                         --run_trio_analysis also mild vs severe promoter methylation per
                                         parental haplotype (<results>/10_modifier_analysis/methylation)
    --run_snpeff [bool]                 snpEff + population AF annotation of the master VCFs (default: true)
    --run_cnv_compare / --run_repeat_compare [bool]   mild vs severe Spectre / Straglr comparison
    --run_read_length_plot [bool]       Read-length histogram + N50 per sample (default: true)

    See nextflow.config for every remaining option (paraphase, 1000 Genomes AF,
    CNV/repeat/methylation thresholds, ...).
    --help                              Show this message
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

params.bin = "${params.base}/../../bin"
params.index = "${params.base}/index"

nextflow.enable.dsl=2

/* modules */
include { DESIGN_INPUT } from '../../modules/parse_design/main.nf'
include { MERGE_INPUTS; SAMTOOLS_CONVERT_BAM_TO_FASTQ; INDEX_REFERENCE } from '../../modules/samtools/main.nf'
include { DORADO_TRIM } from '../../modules/dorado/main.nf'
include { NANOPLOT_RAW } from '../../modules/nanoplot/main.nf'
include { READ_LENGTH_PLOT } from '../../modules/read_length_plot/main.nf'
include { WHATSHAP_PHASE_FAMILY } from '../../modules/whatshap/main.nf'
include { SV_TRIO_ORIGIN } from '../../modules/sv_trio_origin/main.nf'

/* subworkflows */
include { pod5_basecall } from '../../subworkflows/pod5.nf'
include { herro_correction } from '../../subworkflows/herro.nf'
include { trio_family } from '../../subworkflows/trio.nf'
include { assembly } from '../../subworkflows/assembly.nf'
include { minimap2_index; minimap2; minimap2_methylation; clair3; longphase; sniffles; longphase_sv;
          straglr; haplotype_annotate_straglr; spectre; haplotype_annotate_spectre; cnv_compare;
          repeat_compare; paraphase; methylation_haplotag; modkit; consolidate_variants; snpeff;
          kg_af } from '../../subworkflows/ont.nf'
include { haplotype_analysis; modifier_analysis; methylation_compare } from '../../subworkflows/haplotypes.nf'
include { allele_frequency; modifier_annotation } from '../../subworkflows/annotation.nf'

workflow parse_design {
    take :
        data

    main :
        DESIGN_INPUT( data )

        ont = DESIGN_INPUT.out.fastq_ch
            .splitCsv( header: ['sample', 'ont', 'unal_bam'], sep: ",", skip: 1 )
            .map { row -> [ row.sample, file(row.ont), file(row.unal_bam) ] }
            .groupTuple()

        MERGE_INPUTS( ont )

    emit :
        ont = MERGE_INPUTS.out.merged_inputs_ch
}

/*
//////////////////////// MAIN WORKFLOW ////////////////////////
*/

workflow {

    if (!params.pod5_design && !params.design) { exit 1, 'Neither --design nor --pod5_design specified!' }
    if (!params.genome) { exit 1, 'Genome fasta not specified!' }
    ch_genome = file(params.genome, checkIfExists: true)

    run_trio = params.run_trio_analysis && params.family_json
    if (params.run_trio_analysis && !params.family_json) {
        log.warn "--run_trio_analysis is set but --family_json is not -- running without trio analysis"
    }
    if (!params.gnomad_chm13_vcf) {
        log.warn "params.gnomad_chm13_vcf not set -- no gnomAD allele frequencies"
    }
    if (!params.kg_af_vcf && !params.kg_chm13_vcfs) {
        log.info "1000 Genomes CHM13 AF not configured (kg_chm13_vcfs / kg_af_vcf) -- allele frequencies come from gnomAD only"
    }

    // reference indexes: minimap2 .mmi (alignment) and samtools .fai
    minimap2_index( params.genome, params.index )
    reference_index = minimap2_index.out.index
    INDEX_REFERENCE( params.genome )
    samtools_fai_index = INDEX_REFERENCE.out.ref_indexed_ch

    ////////////////// 1. READS: input -> trim -> QC -> HERRO /////////////////////
    if (params.pod5_design) {
        pod5_basecall( params.pod5_design )
        ont_reads = pod5_basecall.out.ont
    } else {
        parse_design( params.design )
        ont_reads = parse_design.out.ont
    }
    if (params.debug) { ont_reads.view() }

    DORADO_TRIM( ont_reads )
    SAMTOOLS_CONVERT_BAM_TO_FASTQ( DORADO_TRIM.out.dorado_trim_output_ch )
    ont_ch = SAMTOOLS_CONVERT_BAM_TO_FASTQ.out.samtool_convert_unalbam2fastq   // (sid, unal_bam, fastq)

    NANOPLOT_RAW( ont_ch )
    if (params.run_read_length_plot) {
        READ_LENGTH_PLOT( ont_ch, file("${params.bin}/plot_read_length.py") )
    }

    // methylation needs the raw (pre-HERRO) reads -- HERRO drops the MM/ML tags
    if (params.run_methylation) {
        minimap2_methylation( reference_index, ont_ch.map { sid, unal_bam, fastq -> tuple(sid, fastq) } )
    }

    herro_correction( ont_ch )
    reads_ch = herro_correction.out.ont_ch

    ////////////////// 2. FAMILY: split parents out /////////////////////
    // Parents only get alignment + Clair3 + Sniffles (inside trio_family);
    // everything below runs on children only.
    if (run_trio) {
        trio_family( params.family_json, reads_ch, reference_index, samtools_fai_index )
        reads_ch = trio_family.out.non_parent_ont_ch
    }

    ////////////////// (optional) ASSEMBLY -- off by default /////////////////////
    // Standalone branch: dipcall / hapdiff calls are published under
    // 04_assembly / 06_variants but are NOT merged into the alignment-based
    // master VCF (stage 6a) -- the two call sets are kept separate.
    if (params.assemble_genome) {
        assembly( reads_ch, run_trio ? trio_family.out.trio_assembly : Channel.empty(),
                  reference_index, samtools_fai_index )
    }

    ////////////////// 3. ALIGNMENT, VARIANT CALLING, PHASING /////////////////////
    minimap2( reference_index, reads_ch.map { sid, unal_bam, fastq -> tuple(sid, fastq) } )
    clair3( samtools_fai_index, minimap2.out.bams )

    // round 1: phase SNVs/indels, haplotag -> Sniffles (phased SVs)
    longphase( samtools_fai_index, clair3.out.clair3_ch.join( minimap2.out.bams ) )
    sniffles( longphase.out.longphase_ch )

    // round 2: re-phase SNVs + SVs together; this BAM/VCF feeds everything below
    longphase_sv( samtools_fai_index,
                  clair3.out.clair3_ch.join( sniffles.out.sniffles_ch ).join( minimap2.out.bams ) )
    phased_bam_ch = longphase_sv.out.longphase_sv_ch          // (sid, bam, bai)
    phased_vcf_ch = longphase_sv.out.longphase_sv_vcf_ch      // (sid, vcf, tbi)

    ////////////////// 4. REPEATS, CNV, PARAPHASE, METHYLATION /////////////////////
    straglr( phased_bam_ch, samtools_fai_index )
    haplotype_annotate_straglr( straglr.out.straglr.join( phased_bam_ch ) )

    spectre( phased_bam_ch, samtools_fai_index )
    haplotype_annotate_spectre( spectre.out.spectre_bed_ch.join( phased_bam_ch ) )

    paraphase( phased_bam_ch, samtools_fai_index )

    if (params.run_methylation) {
        methylation_haplotag( phased_bam_ch, minimap2_methylation.out.bams )
        modkit( methylation_haplotag.out.haplotagged, samtools_fai_index )
    }

    // mild vs severe CNV / repeat comparison
    // (<results>/10_modifier_analysis/{cnv,repeat_expansion})
    repeat_pairs_ch = Channel.value(file("NO_FILE"))
    cnv_gene_cn_ch = Channel.value(file("NO_FILE"))
    if (params.family_json && (params.run_cnv_compare || params.run_repeat_compare)) {
        def cmp_gtf_path = params.cnv_repeat_gtf ?: params.snpeff_gtf
        if (params.run_cnv_compare) {
            if (!cmp_gtf_path) { exit 1, 'params.cnv_repeat_gtf (or params.snpeff_gtf) must be set for run_cnv_compare' }
            cnv_compare( file(params.family_json), spectre.out.spectre_vcf_ch, file(cmp_gtf_path, checkIfExists: true) )
            cnv_gene_cn_ch = cnv_compare.out.cnv_gene_cn_ch
        }
        if (params.run_repeat_compare) {
            repeat_compare(
                file(params.family_json),
                straglr.out.straglr,
                cmp_gtf_path ? file(cmp_gtf_path, checkIfExists: true) : file("NO_FILE"),
                params.repeat_known_loci_bed ? file(params.repeat_known_loci_bed, checkIfExists: true) : file("NO_FILE")
            )
            repeat_pairs_ch = repeat_compare.out.repeat_pairs_ch
        }
    }

    ////////////////// 5. TRIO PHASING + SV PARENT OF ORIGIN /////////////////////
    // small_vcf_ch = the child's phased small-variant VCF used by everything
    // below: WhatsHap pedigree phasing when parents are available, otherwise
    // longphase_sv as-is.
    small_vcf_ch = phased_vcf_ch
    sv_origin_ch = Channel.empty()        // (child, sv_trio_origin.tsv)
    parent_pair_ch = Channel.empty()      // (child, [father_vcf, mother_vcf])
    if (run_trio) {
        parent_pair_ch = trio_family.out.trio_phasing
            .map { child, father, mother, father_vcf, mother_vcf -> tuple(child, [father_vcf, mother_vcf]) }

        // ONE WhatsHap pedigree job per family (all children of a
        // father/mother pair + both parents, with everyone's reads), split
        // back into one phased VCF per child
        parent_bam_map = trio_family.out.parent_bam
            .toList()
            .map { list -> list.collectEntries { pid, bam, bai -> [(pid): bam] } }
        family_phasing_in = phased_vcf_ch.join( phased_bam_ch ).join( trio_family.out.trio_phasing )
            .map { sid, vcf, tbi, bam, bai, father, mother, father_vcf, mother_vcf ->
                tuple("${father}_${mother}".toString(), father, mother, father_vcf, mother_vcf, sid, vcf, bam) }
            .groupTuple( by: 0 )
            .combine( parent_bam_map )
            .map { fam, fathers, mothers, fvcfs, mvcfs, sids, vcfs, bams, pbam ->
                def use_reads = params.whatshap_family_parent_reads
                tuple(fam, fathers[0], mothers[0], fvcfs[0], mvcfs[0],
                      use_reads && pbam[fathers[0]] ? pbam[fathers[0]] : 'NONE',
                      use_reads && pbam[mothers[0]] ? pbam[mothers[0]] : 'NONE',
                      sids, vcfs, bams) }
        WHATSHAP_PHASE_FAMILY( samtools_fai_index, family_phasing_in )
        small_vcf_ch = WHATSHAP_PHASE_FAMILY.out.child_vcfs_ch       // -> (child, vcf, tbi)
            .flatMap { vcfs, tbis ->
                def vl = vcfs instanceof List ? vcfs : [vcfs]
                def tl = tbis instanceof List ? tbis : [tbis]
                vl.collect { v -> tuple(v.name - '.trio.whatshapphase.vcf.gz', v, tl.find { it.name == v.name + '.tbi' }) } }

        // SVs: WhatsHap can't phase them, so their parent of origin comes from
        // Sniffles' own phase tags (assigned on the round-1 longphase BAM, hence
        // round-1's phase blocks) + whether each parent carries the SV.
        // Written into the master VCF's INFO/ORIGIN at consolidation.
        SV_TRIO_ORIGIN(
            sniffles.out.sniffles_ch.map { sid, vcf, tbi -> tuple(sid, vcf) }
                .join( longphase.out.longphase_vcf_ch.map { sid, vcf, tbi -> tuple(sid, vcf) } )
                .join( trio_family.out.trio_sv ),
            file("${params.bin}/sv_trio_origin.py") )
        sv_origin_ch = SV_TRIO_ORIGIN.out.sv_trio_origin_ch.map { sid, tsv, blocks -> tuple(sid, tsv) }
    }

    /*
    ===================== STAGES 6-9 AT A GLANCE =====================
      stage 5 gave us, per child:  phased SNVs/indels (WhatsHap family)
                                   + SVs (Sniffles) + SV parent of origin
      6a CONSOLIDATE  merge all callers -> ONE master VCF per child,
                      every site labelled INFO/ORIGIN (paternal/maternal)
      6b SNPEFF       gene + consequence for every site (ANN)
      6c ALLELE FREQ  gnomAD + 1000 Genomes AF written into INFO, once
      7  IBD          hap-ibd: which parental copies mild & severe sibs share
      8  MILD COPIES  mild variants on a copy the severe sib didn't inherit
      9  METHYLATION  mild vs severe promoter 5mC per parental copy
    ==================================================================
    */

    ////////////////// 6a. CONSOLIDATE CALLERS /////////////////////
    // in:  phased small variants (small_vcf_ch), Sniffles SVs -- alignment-
    //      based callers only (assembly calls stay separate) -- plus parents'
    //      Clair3 VCFs and SV_TRIO_ORIGIN
    // out: <child>.master.vcf -- one row per site, each caller's GT:PS in its
    //      own column, INFO/CALLERS, and (with parents) INFO/ORIGIN
    // why here: WhatsHap only phases SNVs/indels, so this is where they meet
    //      the SVs; ORIGIN needs the trio-phased GTs + SV_TRIO_ORIGIN, both
    //      made in stage 5. Every later step reads this one file per child.
    no_file = file("NO_FILE")
    consolidate_variants(
        small_vcf_ch.map { sid, vcf, tbi -> tuple(sid, vcf) }
            .join( sniffles.out.sniffles_ch.map { sid, vcf, tbi -> tuple(sid, vcf) } )
            .join( parent_pair_ch, remainder: true )
            .join( sv_origin_ch, remainder: true )
            // remainder also emits right-only items (e.g. a family.json child
            // with no reads) -- drop those; index rather than destructure
            .filter { it.size() > 2 && it[1] != null && it[2] != null }
            .map { it ->
                def parents = it.size() > 3 ? it[3] : null
                def sv_tsv  = it.size() > 4 ? it[4] : null
                tuple(it[0], it[1], it[2],
                      parents ? parents[0] : no_file, parents ? parents[1] : no_file, sv_tsv ?: no_file) }
    )

    ////////////////// 6b + 6c. SNPEFF, POPULATION ALLELE FREQUENCY /////////////////////
    // 6b in: master VCF.  out: same VCF + INFO/ANN (gene, consequence, impact).
    // 6c in: snpEff VCFs.  out: <child>.snpeff.af.vcf.gz with INFO/gnomAD_AF_NFE
    //    and INFO/KG_AF_EUR (06_variants/snpeff_af). Two sources because each
    //    misses what the other has: the gnomAD file is a projection onto CHM13
    //    and drops sites where CHM13 carries the minor allele; 1000 Genomes is
    //    re-called natively on CHM13 (smaller, EUR). Annotation only.
    if (params.run_snpeff) {
        if (!params.snpeff_db) { exit 1, 'params.snpeff_db not set!' }
        if (!params.snpeff_data_dir) { exit 1, 'params.snpeff_data_dir not set!' }

        snpeff( consolidate_variants.out.consolidated_ch, ch_genome,
                params.snpeff_gtf ? file(params.snpeff_gtf, checkIfExists: true) : file("NO_FILE") )

        kg_af_vcf_ch = Channel.value(file("NO_FILE"))
        if (params.kg_af_vcf || params.kg_chm13_vcfs) {
            kg_af()
            kg_af_vcf_ch = kg_af.out.kg_af_ch
        }
        allele_frequency( snpeff.out.snpeff_vcf_ch, kg_af_vcf_ch )
    }

    ////////////////// 7. SIBLING IBD WITH HAP-IBD /////////////////////
    // in:  WHATSHAP_PHASE_FAMILY's family VCF (all members, one per family)
    // ORIENT_PHASE: children only, each child's phase blocks flipped so GT is
    //      paternal|maternal along the whole chromosome (hap-ibd ignores PS)
    // HAP_IBD: haplotype-level IBD segments (uniform 1 cM/Mb map)
    // IBD_SEGMENTS: per mild x severe pair, per parent, shared / not_shared,
    //      combined into IBD0/1/2 (<results>/07_haplotypes)
    // Runs on the family VCF, not the master VCFs: hap-ibd needs every
    // member's phased genotypes in one VCF and only uses biallelic SNVs.
    if (run_trio) {
        haplotype_analysis( file(params.family_json), WHATSHAP_PHASE_FAMILY.out.family_vcf_ch, samtools_fai_index )

        ////////////////// 8. MODIFIER GENES (variant based) /////////////////////
        // in:  annotated master VCFs (ORIGIN + ANN + AF) + the stage 7 IBD0/1/2
        //      segments + the gene model
        // out: per family (<results>/10_modifier_analysis/variants/families/<family>):
        //      the mild child's passing genes -- homozygous (severe sib not
        //      homozygous for it), compound_het (>= 1 het per parental copy,
        //      >= 1 absent from the severe sib), regulatory (non-genic, not
        //      carried identically by the severe sib; both flanking genes) --
        //      plus mild / severe VCFs of the supporting variants with
        //      INFO/IBD + AF. Then the cohort count of families per gene.
        //      No statistical model and no AF filter.
        if (params.run_modifier_analysis) {
            if (!params.run_snpeff) {
                log.warn "--run_modifier_analysis needs snpEff annotation (--run_snpeff) -- skipping the modifier search"
            } else {
                def mod_gtf_path = params.cnv_repeat_gtf ?: params.snpeff_gtf
                if (!mod_gtf_path) { exit 1, 'params.cnv_repeat_gtf (or params.snpeff_gtf) must be set for run_modifier_analysis' }
                modifier_analysis(
                    file(params.family_json), allele_frequency.out.af_vcf,
                    haplotype_analysis.out.ibd_segments_ch.first(),
                    file(mod_gtf_path, checkIfExists: true)
                )

                ////////////////// 10. ALPHAGENOME / ALPHAMISSENSE /////////////////////
                // every mild / severe variant in the genes passing in >=
                // modifier_min_mild mild children, one chain per side; each
                // tool's per-variant output as is
                // (<results>/10_modifier_analysis/variants/<tool>/{mild,severe})
                if (params.run_alphagenome || params.run_alphagenome_api || params.alphamissense_tsv) {
                    modifier_annotation(
                        modifier_analysis.out.mild_candidate_vcf_ch,
                        modifier_analysis.out.severe_candidate_vcf_ch,
                        samtools_fai_index
                    )
                }
            }
        }

        ////////////////// 9. PROMOTER METHYLATION /////////////////////
        // in:  modkit HP1/HP2 pileups, longphase_sv VCF (source of the HP
        //      tags), parents' VCFs, the stage 7 sharing segments
        // out: mild vs severe promoter 5mC, combined and per parental copy,
        //      each copy labelled shared / not_shared (<results>/10_modifier_analysis/methylation)
        if (params.run_methylation) {
            def meth_gtf_path = params.cnv_repeat_gtf ?: params.snpeff_gtf
            if (!meth_gtf_path) { exit 1, 'params.cnv_repeat_gtf (or params.snpeff_gtf) must be set for run_methylation' }
            methylation_compare(
                file(params.family_json), modkit.out.modkit, phased_vcf_ch,
                trio_family.out.parent_vcf, haplotype_analysis.out.parent_segments_ch,
                file(meth_gtf_path, checkIfExists: true)
            )
        }
    }
}
