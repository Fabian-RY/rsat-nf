process DYAD{

    container "biocontainers/rsat:2025-03-26_cv1"

    input:
    tuple val(meta), path(input)
    val organism

    output:
    tuple val(meta), path("*.txt")       

    script:
    """
        dyad-analysis -i ${input} -v 1 -quick -sort -timeout 3600  -type any -2str -noov -lth occ 1 -lth occ_sig 0 -uth rank 50 -return occ,proba,rank -l 3 -spacing 0-20  -bg upstream-noorf -org ${organism} -o ${meta}.txt
    """

}