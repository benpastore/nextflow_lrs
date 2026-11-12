nextflow.enable.dsl=2

// De novo assembly branch (hifiasm -> gfatools -> dorado polish -> dipcall /
// hapdiff). Moved out of main.nf unchanged; OFF by default
// (params.assemble_genome = false) while the analysis focuses on the
// alignment-based calls. Turn it back on with --assemble_genome true and its
// dipcall/hapdiff calls rejoin consolidate_variants automatically.

include { hifasim; hifasim_trio; gfatools; dorado; genome_stats; minimap2_map_asm_to_ref;
          chrom_coverage; dipcall; hapdiff } from './ont.nf'

workflow assembly {

    take :
        ont_ch               // tuple(sampleID, unal_bam, fastq) -- herro-corrected, children only
        trio_assembly        // tuple(child, paternal_yak, maternal_yak), or Channel.empty() without trios
        reference_index
        samtools_fai_index

    main :
        // children with both parents' yak dbs get trio-binning, everyone else
        // the regular assembly. join(remainder) pads a miss with ONE null, so
        // no-trio tuples are (sid, bam, fastq, null) -- index, don't destructure.
        gate = ont_ch
            .join( trio_assembly, remainder: true )
            .branch {
                has_trio: it.size() > 3 && it[3] != null
                no_trio: true
            }
        hifasim( gate.no_trio.map { it -> tuple(it[0], it[1], it[2]) } )
        hifasim_trio( gate.has_trio.map { sid, bam, fastq, pat_yak, mat_yak -> tuple(sid, bam, fastq, pat_yak, mat_yak) } )
        asm_ch = hifasim.out.hifasim_asm.mix( hifasim_trio.out.hifasim_trio_asm )

        gfatools( asm_ch )

        if (params.dorado) {
            dorado( gfatools.out.fasta_asm )
            genome_asm_ch = dorado.out.dorado_output_ch
        } else {
            genome_asm_ch = gfatools.out.haplotype_asm
                .map { sampleID, bam, reads, hap1, hap2 -> tuple(sampleID, reads, hap1, hap2) }
        }

        genome_stats( genome_asm_ch )
        minimap2_map_asm_to_ref( reference_index, genome_asm_ch )
        chrom_coverage( minimap2_map_asm_to_ref.out.bams )

        dipcall( genome_asm_ch, samtools_fai_index )
        hapdiff( genome_asm_ch, samtools_fai_index )

    emit :
        dipcall_vcf = dipcall.out.dipcall.map { sampleID, vcf, tbi -> tuple(sampleID, vcf) }
        hapdiff_vcf = hapdiff.out.hapdiff_ch
}
