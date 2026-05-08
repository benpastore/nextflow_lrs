/*
========================================================================================
                            Verkko Assembly Module
========================================================================================
Verkko v1.4.1 - Telomere-to-telomere genome assembler
https://github.com/marbl/verkko
*/

process VERKKO {
    
    tag "${sample}"
    label 'verkko'
    
    publishDir "${params.results}/assembly/verkko/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), path(hifi_reads), path(ont_reads)
    
    output:
    tuple val(sample), path("${sample}.verkko.fasta"), emit: assembly
    path "${sample}.verkko.gfa", emit: graph
    path "${sample}_verkko_stats.txt", emit: stats
    path "verkko_output/", emit: full_output
    
    script:
    def ont_arg = ont_reads ? "--nano ${ont_reads}" : ""
    """
    # Run Verkko
    verkko \\
        -d verkko_output \\
        --hifi ${hifi_reads} \\
        ${ont_arg} \\
        --threads ${task.cpus} \\
        ${params.verkko_extra_args}
    
    # Extract assembly
    cp verkko_output/assembly.fasta ${sample}.verkko.fasta
    
    # Extract graph if available
    if [ -f verkko_output/assembly.gfa ]; then
        cp verkko_output/assembly.gfa ${sample}.verkko.gfa
    else
        touch ${sample}.verkko.gfa
    fi
    
    # Generate stats
    echo "Sample: ${sample}" > ${sample}_verkko_stats.txt
    echo "Assembly Statistics:" >> ${sample}_verkko_stats.txt
    
    # Get contig count and N50
    if command -v seqkit &> /dev/null; then
        seqkit stats ${sample}.verkko.fasta >> ${sample}_verkko_stats.txt
    else
        # Simple stats without seqkit
        num_contigs=\$(grep -c "^>" ${sample}.verkko.fasta)
        echo "Number of contigs: \$num_contigs" >> ${sample}_verkko_stats.txt
    fi
    """
}