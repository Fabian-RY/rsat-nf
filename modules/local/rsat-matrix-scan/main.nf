process MATRIX_SCAN {
    
    container "biocontainers/rsat:2025-03-26_cv1"

    input:
    tuple val(meta), path(input)
    val input

    script:
    """
    matrix-scan -v 1 -matrix_format transfac -m ${SCAN_INPUT}/oligos_5-8nt_m1/peaks_oligos_5-8nt_m1.tf -i ${UP1}/regulon${REGULON}_up1.rm.fna \
    -seq_format fasta \
    -pseudo 1 -decimals 1 -2str -origin end -bginput -markov 1 -bg_pseudo 0.01 -return limits -return sites -return pval -lth score 1 -uth pval ${SCANMAXPVALUE} \
    -n score -o ${SCAN_OLIGO_OUTPUT}/scan_oligo_up1.tabecho "test matrix_scan"
    """

}
