process PEAK_MOTIFS{
    
    container "biocontainers/rsat:2025-03-26_cv1"

    input:
    tuple val(from_to), val(meta), path(motifs), val(meta2), path(control)
    val motif_db
    val footdbfile

    output:
    tuple val(meta), val(from_to), emit: outdir

    script:
    """
         peak-motifs -i ${motifs} -ctrl ${control} -motif_db ${motif_db} transfac ${footdbfile} -prefix peaks -outdir ${meta}_${from_to}_outdir -title analysis_M11 -origin end -disco oligos,dyads -nmotifs 5 -minol 5 -maxol 8 -scan_markov 1 -noov -img_format png -task purge,seqlen,composition,disco,merge_motifs,split_motifs,timelog,motifs_vs_db,synthesis,small_summary,clean_seq
    """
}

process PEAK_MOTIFS_RANDOM{
    
    container "biocontainers/rsat:2025-03-26_cv1"

    input:
    tuple val(from_to), val(meta), path(motifs), val(meta2), path(control)
    val motif_db
    val footdbfile

    output:
    tuple val(meta), path("*.pk"), emit: outdir

    script:
    """
         peak-motifs -i ${motifs} -ctrl ${control} -motif_db ${motif_db} transfac ${footdbfile} -prefix peaks -outdir . -title analysis_M11 -origin end -disco oligos,dyads -nmotifs 5 -minol 5 -maxol 8 -scan_markov 1 -noov -img_format png -task purge,seqlen,composition,disco,merge_motifs,split_motifs,timelog,motifs_vs_db,synthesis,small_summary,clean_seq
    """
}