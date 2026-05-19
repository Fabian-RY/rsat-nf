process RANDOM_GENES {

    container "biocontainers/rsat:2025-03-26_cv1"

    input:
    val number
    val organism
    val feature_type

    output:
    tuple val(number), path("out.txt"), emit: outfile

    script:
    """
        random-genes -n ${number} -org ${organism} -feattype ${feature_type} -o out.txt
    """
}