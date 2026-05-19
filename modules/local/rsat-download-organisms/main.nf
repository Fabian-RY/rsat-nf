process DOWNLOAD_ORGANISM {
    container "biocontainers/rsat:2025-03-26_cv1"

    input:
    val server
    val organism

    output:
    tuple val(organism), val(server), emit: organism

    script:
    """
        download-organism -v 2 -server ${server} -org ${organism}
    """

}