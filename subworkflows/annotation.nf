nextflow.enable.dsl=2

// Annotation of modifier candidates + final report. Moved out of main.nf
// unchanged: takes the cohort recessive/dominant TSVs (from either the
// haplotype-aware or the legacy modifier search -- same schema) and adds
// AlphaGenome, AlphaMissense and the per-variant report.

include { INDEX_REFERENCE as INDEX_HG38_REFERENCE } from '../modules/samtools/main.nf'
include { alphagenome; alphamissense; variant_report } from './ont.nf'

workflow candidate_annotation {

    take :
        recessive_tsv
        dominant_tsv
        family_json
        vcf_manifest
        vcfs
        samtools_fai_index
        repeat_pairs        // or file("NO_FILE")
        cnv_gene_cn         // or file("NO_FILE")

    main :
        hg38_vcf_ch = Channel.empty()
        candidate_vcf_ch = Channel.empty()
        unmap_ch = Channel.empty()
        alphagenome_ch = Channel.empty()
        alphagenome_api_ch = Channel.empty()
        alphagenome_gene_ch = Channel.empty()
        alphagenome_api_gene_ch = Channel.empty()
        alphamissense_ch = Channel.empty()

        // AlphaGenome / AlphaMissense are hg38 tools: lift CHM13 candidates
        // over (genome_version T2T) or use params.genome as-is (HG38)
        if (params.run_alphagenome || params.run_alphagenome_api || params.alphamissense_tsv) {
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
        }

        if (params.run_alphagenome || params.run_alphagenome_api) {
            if (params.run_alphagenome && !params.alphagenome_weights) { exit 1, 'params.alphagenome_weights not set!' }
            alphagenome( recessive_tsv, dominant_tsv, hg38_fai_index, skip_liftover )
            hg38_vcf_ch = alphagenome.out.hg38_vcf_ch
            candidate_vcf_ch = alphagenome.out.candidate_vcf_ch
            unmap_ch = alphagenome.out.alphagenome_unmapped_ch
            alphagenome_ch = alphagenome.out.alphagenome_ch
            alphagenome_api_ch = alphagenome.out.alphagenome_api_ch
            alphagenome_gene_ch = alphagenome.out.alphagenome_gene_ch
            alphagenome_api_gene_ch = alphagenome.out.alphagenome_api_gene_ch
        }

        if (params.alphamissense_tsv) {
            alphamissense( recessive_tsv, dominant_tsv, hg38_fai_index, skip_liftover )
            alphamissense_ch = alphamissense.out.alphamissense_ch
        }

        // one row per candidate; whatever wasn't run is left as empty columns
        variant_report(
            recessive_tsv, dominant_tsv, family_json, vcf_manifest, vcfs,
            candidate_vcf_ch.ifEmpty(file("NO_FILE")),
            hg38_vcf_ch.ifEmpty(file("NO_FILE")),
            unmap_ch.ifEmpty(file("NO_FILE")),
            alphagenome_ch.ifEmpty(file("NO_FILE")),
            alphagenome_api_ch.ifEmpty(file("NO_FILE")),
            alphamissense_ch.ifEmpty(file("NO_FILE")),
            repeat_pairs,
            cnv_gene_cn,
            alphagenome_gene_ch.ifEmpty(file("NO_FILE")),
            alphagenome_api_gene_ch.ifEmpty(file("NO_FILE"))
        )
}
