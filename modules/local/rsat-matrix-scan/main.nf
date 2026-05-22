process MATRIX_SCAN {
    
    container "biocontainers/rsat:2025-03-26_cv1"

    input:
    tuple val(meta), path(input)
    val organism

    script:
    """
    echo "test matrix_scan"
    """

}
