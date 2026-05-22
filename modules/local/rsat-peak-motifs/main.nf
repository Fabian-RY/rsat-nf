process PEAK_MOTIFS{
    
    container "biocontainers/rsat:2025-03-26_cv1"

    input:
    tuple val(meta),  path(motifs)
    tuple val(meta2), path(control)

    output:
    tuple val(meta), path("*.pk"), emit: outdir

    script:
    """
        peak-motifs -v 1 -title ${meta} -i ${motifs} -ctrl ${control} -markov auto -disco oligos -nmotifs 5  -minol 6 -maxol 7  -no_merge_lengths     -2str      -origin center      -scan_markov 1     -source galaxy      -task purge,seqlen,composition,disco,merge_motifs,split_motifs,motifs_vs_motifs,timelog,archive,synthesis,small_summary,scan     -prefix peak-motifs     -noov     -img_format png      -outdir ${meta}.pk    
    """
}