nextflow.enable.dsl=2

// allele_frequency: population AF into every child's snpEff VCF.
//
// modifier_annotation: AlphaGenome + AlphaMissense over every mild and
// severe variant in the recurrent modifier genes (bin/gene_recurrence.py's
// candidate_variants.{mild,severe}.chm13.vcf), each side scored by its own
// chain and published under .../variants/<tool>/<side>. INFO/MILD_SAMPLES,
// SEVERE_SAMPLES and CARRIERS are carried into every output row.

include { INDEX_REFERENCE as INDEX_HG38_REFERENCE } from '../modules/samtools/main.nf'
include { alphagenome as alphagenome_mild; alphagenome as alphagenome_severe;
          alphamissense as alphamissense_mild; alphamissense as alphamissense_severe } from './ont.nf'
include { ANNOTATE_AF } from '../modules/af_annotate/main.nf'

workflow modifier_annotation {

    take :
        mild_vcf             // candidate_variants.mild.chm13.vcf
        severe_vcf           // candidate_variants.severe.chm13.vcf
        samtools_fai_index

    main :
        alphagenome_ch = Channel.empty()
        alphagenome_api_ch = Channel.empty()
        alphamissense_ch = Channel.empty()

        // AlphaGenome / AlphaMissense are hg38 tools: lift CHM13 candidates
        // over (genome_version T2T) or use params.genome as-is (HG38)
        if (!(params.genome_version in ['T2T', 'HG38'])) {
            exit 1, "params.genome_version must be 'T2T' or 'HG38' (got: ${params.genome_version})"
        }
        skip_liftover = (params.genome_version == 'HG38')
        if (skip_liftover) {
            hg38_fai_index = samtools_fai_index
        } else {
            if (!params.hg38_genome) { exit 1, 'params.hg38_genome not set! (required when params.genome_version is "T2T")' }
            INDEX_HG38_REFERENCE( file(params.hg38_genome, checkIfExists: true) )
            hg38_fai_index = INDEX_HG38_REFERENCE.out.ref_indexed_ch
        }

        if (params.run_alphagenome || params.run_alphagenome_api) {
            if (params.run_alphagenome && !params.alphagenome_weights) { exit 1, 'params.alphagenome_weights not set!' }
            alphagenome_mild( 'mild', mild_vcf, hg38_fai_index, skip_liftover )
            alphagenome_severe( 'severe', severe_vcf, hg38_fai_index, skip_liftover )
            alphagenome_ch = alphagenome_mild.out.alphagenome_ch.mix( alphagenome_severe.out.alphagenome_ch )
            alphagenome_api_ch = alphagenome_mild.out.alphagenome_api_ch.mix( alphagenome_severe.out.alphagenome_api_ch )
        }

        if (params.alphamissense_tsv) {
            alphamissense_mild( 'mild', mild_vcf, hg38_fai_index, skip_liftover )
            alphamissense_severe( 'severe', severe_vcf, hg38_fai_index, skip_liftover )
            alphamissense_ch = alphamissense_mild.out.alphamissense_ch.mix( alphamissense_severe.out.alphamissense_ch )
        }

    emit :
        alphagenome_ch     = alphagenome_ch
        alphagenome_api_ch = alphagenome_api_ch
        alphamissense_ch   = alphamissense_ch
}

// Population AF written into every child's snpEff VCF once, right after
// snpEff (see modules/af_annotate/main.nf).

workflow allele_frequency {

    take :
        snpeff_vcf   // tuple(sampleID, snpEff VCF)
        kg_af_vcf    // 1000 Genomes CHM13 EUR sites VCF, or file("NO_FILE")

    main :
        manifest = snpeff_vcf
            .map { sid, vcf -> "${sid}\t${vcf.name}" }
            .collectFile(name: 'af_vcf_manifest.tsv', newLine: true, sort: true)
        ANNOTATE_AF(
            manifest, snpeff_vcf.map { sid, vcf -> vcf }.collect(), kg_af_vcf,
            file("${params.bin}/_modifier_common.py"),
            file("${params.bin}/annotate_af.py")
        )
        af_vcf = ANNOTATE_AF.out.af_vcfs_ch
            .flatMap { vcfs -> (vcfs instanceof List ? vcfs : [vcfs])
                .collect { v -> tuple(v.name - '.snpeff.af.vcf.gz', v) } }

    emit :
        af_vcf = af_vcf   // tuple(sampleID, <sampleID>.snpeff.af.vcf.gz)
}
