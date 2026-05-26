process RETRIEVE_SEQUENCES {

    container "biocontainers/rsat:2025-03-26_cv1"
    
    input:
    val organism
    tuple val(meta), path(module), 
    val(meta_ch), val(from), val(to)
    val output
    val feattype
    val type
    val format
    val label

    output:
    tuple val(meta), path("*.fasta"), val("${from}_${to}")

    script:
    """
    retrieve-seq -org ${organism} -feattype ${feattype} -format ${format} -label ${label} -from ${from} -to ${to} -type ${type} -o ${meta}.fasta -i ${module}
    """

}

process RETRIEVE_GENOME_SEQUENCES {

    container "biocontainers/rsat:2025-03-26_cv1"
    
    input:
    val organism
    tuple val(meta), val(from), val(to)
    val output
    val feattype
    val type
    val format
    val label
    //val gene_list

    output:
    tuple val(meta), path("*.fasta"), val("${from}_${to}"), emit: seqs

    script:
    """
    retrieve-seq -org ${organism} -feattype ${feattype} -from ${from} -to ${to} -noorf -all -label id -rm -o ${organism}_${from}_${to}.fasta
    """

}

process REPEAT_CHANNEL_SEQUENCES {

    input:
    tuple val(meta), path(bg), val(number)

    output:
    tuple val(meta), path(repeated)

    exec:
    repeated = bg * number

}