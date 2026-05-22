process RANDOM_GENES {

    container "biocontainers/rsat:2025-03-26_cv1"

    input:
    val id
    val number
    val organism
    val feature_type

    output:
    tuple val(id), path("*.txt"), emit: outfile

    script:
    """
        random-genes -n ${number} -org ${organism} -feattype ${feature_type} -o ${id}.txt
    """
}