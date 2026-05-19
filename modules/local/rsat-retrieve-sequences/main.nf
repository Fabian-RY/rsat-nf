process RETRIEVE_SEQUENCES {

    container "biocontainers/rsat:2025-03-26_cv1"
    
    input:
    val organism

    output:
    tuple val(organism), path("*.fasta")

    script:
    """

    """

}
