process OLIGO {
    
    container "biocontainers/rsat:2025-03-26_cv1"

    input:
    tuple val(meta), path(input)
    val organism

    output:
    tuple val(meta), path("*.tab")   

    script:
    """
    oligo-analysis  -v 1 -sort -i ${input} -format fasta  -lth occ_sig 0 -uth rank 50 -return occ,proba,rank -2str -noov -quick_if_possible  -seqtype dna -bg upstream -org ${organism} -pseudo 0.01 -l 6 -o ${meta}.tab; 
    """

}