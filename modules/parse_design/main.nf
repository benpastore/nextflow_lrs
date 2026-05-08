
process DESIGN_INPUT {

    label 'low'

    publishDir "$params.results/samples", mode : 'copy', pattern : "*csv"

    input : 
        val(design)

    output : 
        path("fastq.csv"), emit : fastq_ch
        //path("replicates.csv"), emit : condition_ch

    script : 
    """
    #!/bin/bash
        
    source activate rnaseq
    
    python3 ${params.bin}/parse_design.py -input ${design}
    """

}