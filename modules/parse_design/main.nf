
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

    # conda's own activate.d hooks for this env (specifically
    # binutils_linux-64's, which references $ADDR2LINE with no default)
    # aren't written to survive `set -u` -- drop it just for activation.
    set +u
    source activate rnaseq
    set -u

    python3 ${params.bin}/parse_design.py -input ${design}
    """

}