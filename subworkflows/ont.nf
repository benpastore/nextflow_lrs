nextflow.enable.dsl=2

include { HIFASIM } from '../modules/hifasim/main.nf'
workflow hifasim { 

    take :
        data

    main : 
        HIFASIM( data )

    emit : 
        hifasim_asm = HIFASIM.out.hifasim_asm

}

include { GFA_CONVERT } from '../modules/gfatools/main.nf'
workflow gfatools { 
    
    take : 
        data 
    
    main : 
        GFA_CONVERT( data )
    
    emit : 
        fasta_asm = GFA_CONVERT.out.fasta_asm
}

include { MAP_READS_TO_ASSEMBLY } from '../modules/minimap2/main.nf'
workflow map_reads_to_assembly {

    take : 
        data 
    
    main : 
        MAP_READS_TO_ASSEMBLY( data )
    
    emit : 
        assembly_alignment = MAP_READS_TO_ASSEMBLY.out.alignment

}

include { HAPDUP_PHASE } from '../modules/hapdup/main.nf'
workflow hapdup_phase { 

    take : 
        data 
    
    main : 
        HAPDUP_PHASE( data )
    
    emit : 
        hapdup_output = HAPDUP_PHASE.out.hapdup_ch

}


include { DIPCALL } from '../modules/dipcall/main.nf'
workflow dipcall { 

    take : 
        data 
        genome
    
    main : 
        DIPCALL( data, genome )
    
    emit : 
        dipcall = DIPCALL.out.dipcall_ch

}

include { MINIMAP2_INDEX } from '../modules/minimap2/main.nf'
include { MINIMAP2_ALIGN } from '../modules/minimap2/main.nf'
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

include { CLAIR3 } from '../modules/clair3/main.nf'
workflow clair3 { 
    take : 
        data 
    main : 
        CLAIR3( params.genome, data )
    emit : 
        clair3_ch = CLAIR3.out.clair3_ch
}

include { WHATSHAP_PHASE } from '../modules/whatshap/main.nf'
include { WHATSHAP_HAPLOTAG } from '../modules/whatshap/main.nf'
workflow whatshap_phase { 
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

include { SNIFFLES } from '../modules/sniffles/main.nf'
workflow sniffles { 

    take : 
        data 
    main : 
        SNIFFLES( data )
    emit : 
        sniffles = SNIFFLES.out.sniffles_ch

}

include { SPECTRE } from '../modules/spectre/main.nf'
workflow spectre { 

    take : 
        data 
    
    main : 
        SPECTRE( data )
    
    emit : 
        spectre = SPECTRE.out.spectre_ch

}

include { STRAGLR } from '../modules/straglr/main.nf'
workflow straglr { 

    take : 
        data 
    
    main : 
        STRAGLR( data )
    
    emit : 
        straglr = STRAGLR.out.straglr_ch



}