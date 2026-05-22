process RETRIEVE_SEQUENCES {

    container "biocontainers/rsat:2025-03-26_cv1"
    
    input:
    val organism
    tuple val(meta), path(module)
    val output
    val feattype
    val type
    val format
    val label
    val from
    val to
    //val gene_list

    output:
    tuple val(meta), path("*.fasta"), emit: seqs

    script:
    """
    retrieve-seq -org ${organism} -feattype ${feattype} -format ${format} -label ${label} -from ${from} -to ${to} -type ${type} -o ${meta}.fasta -i ${module}
    """

}

process RETRIEVE_ALL_SEQUENCES {

    container "biocontainers/rsat:2025-03-26_cv1"
    
    input:
    val organism
    val output
    val feattype
    val type
    val format
    val label
    val from
    val to
    //val gene_list

    output:
    tuple val(organism), path(output), emit: seqs

    script:
    """
    retrieve-seq -org ${organism} -feattype ${feattype} -format ${format} -label ${label} -from ${from} -to ${to} -type ${type} -o ${output} -all
    """

}
