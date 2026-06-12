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
include { GFA_FAIDX } from '../modules/gfatools/main.nf'
workflow gfatools { 
    
    take : 
        data 
    
    main : 
        GFA_CONVERT( data )
        GFA_FAIDX( GFA_CONVERT.out.combined_fa )
    
    emit : 
        fasta_asm = GFA_FAIDX.out.fasta_asm
        haplotype_asm = GFA_FAIDX.out.haplotype_fasta_asm
}

include { DORADO } from '../modules/dorado/main.nf'
workflow dorado { 

    take : 
        data 
    
    main : 
        DORADO( data )
    
    emit : 
        dorado_output_ch = DORADO.out.dorado_output_ch

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

include { HAPDIFF } from '../modules/hapdiff/main.nf'
workflow hapdiff { 

    take : 
        data 
        genome
    
    main : 
        HAPDIFF( data, genome )
    
    emit : 
        hapdiff_ch = HAPDIFF.out.hapdiff_ch

}

include { MINIMAP2_INDEX } from '../modules/minimap2/main.nf'
include { MINIMAP2_ALIGN } from '../modules/minimap2/main.nf'
workflow minimap2_index { 

    take : 
        data
        index_dir
    
    main : 

        genome_fasta = file("${data}")
        genome_name = "${data.baseName}"
        index_dir = "${index_dir}/minimap2"
        index_path = "${index_dir}/minimap2/${genome_name}.mmi"

        index_exists = file(index_path).exists()

        if ( index_exists == false || params.force_build_mmi){
            build = true
        } else {
            build = false
        }

        if (build == true){
            MINIMAP2_INDEX( data, index_dir, index_path )
            index = MINIMAP2_INDEX.out.index_ch
        } else {
            index = index_path 
        }
    
    emit : 
        index = index
    
}

workflow minimap2 {

    take : 
        index
        data

    main :

        MINIMAP2_ALIGN( index, data )
    
    emit :
        bams = MINIMAP2_ALIGN.out.bam_ch
}

include { MAP_ASM_TO_REF } from '../../modules/minimap2/main.nf'
workflow minimap2_map_asm_to_ref {

    take : 
        index 
        data
        
    main :
        MAP_ASM_TO_REF( index, data )
    
    emit :
        bams = MINIMAP2_ALIGN.out.bam_ch
}


include { CLAIR3 } from '../modules/clair3/main.nf'
workflow clair3 { 
    take : 
        ref 
        data

    main : 
        CLAIR3( ref, data )
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

include { LONGPHASE } from '../modules/longphase/main.nf'
workflow longphase { 

    take : 
        ref
        data
    main : 
        LONGPHASE( ref, data )
    emit : 
        longphase_ch = LONGPHASE.out.longphase_ch        


}

include { LONGPHASE_SV } from '../modules/longphase/main.nf'
workflow longphase_sv { 

    take : 
        ref
        data

    main : 
        LONGPHASE_SV( ref, data )

    emit : 
        longphase_sv_ch = LONGPHASE_SV.out.longphase_sv_ch        


}

include { SNIFFLES } from '../modules/sniffles/main.nf'
include { INDEX_SNIFFLES_VCF } from '../modules/sniffles/main.nf'
workflow sniffles { 

    take : 
        data
        
    main : 
        SNIFFLES( data )

    emit : 
        sniffles_ch = SNIFFLES.out.sniffles_ch

}

include { SPECTRE } from '../modules/spectre/main.nf'
workflow spectre { 

    take : 
        data 
        genome
    
    main : 
        SPECTRE( data, genome )
    
    emit : 
        spectre = SPECTRE.out.spectre_cnv_ch

}

include { STRAGLR } from '../modules/straglr/main.nf'
workflow straglr { 

    take : 
        data 
        genome
    
    main : 
        STRAGLR( data, genome )
    
    emit : 
        straglr = STRAGLR.out.straglr_ch



}