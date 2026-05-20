process PURGE_SEQUENCES {

    container "biocontainers/rsat:2025-03-26_cv1"


    input:
    tuple val(meta), path(input)
    val format

    output:
    tuple val(meta), path("*.purged"), emit: purged

    script:
    """
        purge-sequence -i ${input} -format ${format} -o ${input}.purged 
    """

}