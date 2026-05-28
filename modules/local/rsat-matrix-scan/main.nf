process MATRIX_SCAN {
    
    container "biocontainers/rsat:2025-03-26_cv1"

    input:
    tuple val(from_to), val(meta), path(input), val(up) , val(regulon_fna)
    val max_p_val

    output:
    path("*.tabecho")

    script:
    """
    matrix-scan -v 1 -matrix_format transfac -m ${input} -i ${regulon_fna} -seq_format fasta -pseudo 1 -decimals 1 -2str -origin end -bginput -markov 1 -bg_pseudo 0.01 -return limits -return sites -return pval -lth score 1 -uth pval ${max_p_val} -n score -o ./${meta}_scan.tab
    """

}
