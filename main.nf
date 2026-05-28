#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW_READS} from './modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_TRIMMERED_READS  } from './modules/nf-core/fastqc/main.nf'
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main.nf'
include { KALLISTO_INDEX } from './modules/nf-core/kallisto/index/main.nf' 
include { KALLISTO_QUANT } from './modules/nf-core/kallisto/quant/main.nf'
include { MULTIQC } from './modules/nf-core/multiqc/main.nf'
include { DOWNLOAD_ORGANISM } from './modules/local/rsat-download-organisms/main.nf'
include { RANDOM_GENES } from './modules/local/rsat-random-genes/main.nf'
include { RETRIEVE_SEQUENCES } from './modules/local/rsat-retrieve-sequences/main.nf'
include { RETRIEVE_SEQUENCES as RETRIEVE_RANDOM_SEQUENCES } from './modules/local/rsat-retrieve-sequences/main.nf'
include { RETRIEVE_GENOME_SEQUENCES } from './modules/local/rsat-retrieve-sequences/main.nf'
include { REPEAT_CHANNEL_SEQUENCES } from './modules/local/rsat-retrieve-sequences/main.nf'
include { PURGE_SEQUENCES } from './modules/local/rsat-purge-sequences/main.nf'
include { DYAD } from './modules/local/rsat-dyad-analysis/main.nf'
include { OLIGO } from './modules/local/rsat-oligo-analysis/main.nf'
include { PEAK_MOTIFS} from './modules/local/rsat-peak-motifs/main.nf'
include { PEAK_MOTIFS as PEAK_MOTIFS_RANDOM } from './modules/local/rsat-peak-motifs/main.nf'
include { MATRIX_SCAN } from './modules/local/rsat-matrix-scan/main.nf'
include { MATRIX_SCAN as MATRIX_SCAN_RANDOM } from './modules/local/rsat-matrix-scan/main.nf'

// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

    
    modules_ch = Channel.fromPath(params.modules)
        .ifEmpty { "Error: No samples found in: ${params.samplesheet}"}
        .splitCsv(header: true)
        .map { row -> 
            def meta = [id: row.id.trim()]
            def module = file(row.module.trim())
            [meta, module]    
        }

    background_ch = Channel.fromPath(params.background_sheet)
        .ifEmpty { "Error: No samples found in: ${params.background_sheet}"}
        .splitCsv(header: true)
        .map { row -> 
            def meta = [id: row.id.trim()]
            def from = row.from.trim()
            def to = row.to.trim()
            [meta, from, to]    
        }


    all_genes = Channel.of(tuple(params.organism, "-all"))

    // Checks to download the organism or not
    if (params.download_organism){
        //download_result = DOWNLOAD_ORGANISM(params.server, params.organism)
        //def organism_downloaded = download_result.map {org, server -> org}
        //genes = RANDOM_GENES("random", params.number, organism_downloaded, params.feature_type)
        //sequences = RETRIEVE_ALL_SEQUENCES(organism_downloaded, params.retrieve_seq_output,
        //        params.feature_type, params.retrieve_seq_type, params.retrieve_seq_format, params.retrieve_seq_label, params.retrieve_seq_from,
        //        params.retrieve_seq_to)
    }
    else {

        ///////////////
        //
        // First part of pipeline: Retrieving sequences
        //
        // We need to retrieve: 
        // - Background sequences 
        // - Motif Sequences 
        // - Random sequences


        // Retrieve all background sequences
        background_sequences = RETRIEVE_GENOME_SEQUENCES(params.organism, background_ch,
                params.retrieve_seq_type, params.feature_type, params.retrieve_seq_format, params.retrieve_seq_label, background_ch)
        
        
        combinations_ch = modules_ch.combine(background_ch)
        // Retrieve module sequences in this boundaries
        regulon_sequences = RETRIEVE_SEQUENCES(params.organism, combinations_ch, params.retrieve_seq_output,
                                    params.feature_type, params.retrieve_seq_type, params.retrieve_seq_format, params.retrieve_seq_label)
        // Generate random modules and Retrieve random modules
        n_clusters_ch = Channel
                        .from(1..params.clusters)
                        .each {
                            i -> "random_cluster_${1}"
                        }
        clusters = RANDOM_GENES(n_clusters_ch, params.number, params.organism, params.feature_type)         
        random_combinations = clusters.combine(background_ch)
        random_sequences = RETRIEVE_RANDOM_SEQUENCES(params.organism, random_combinations, params.retrieve_seq_output,
                params.feature_type, params.retrieve_seq_type, params.retrieve_seq_format, params.retrieve_seq_label)

        /////////////////////////////////////////////
        //
        // PEAK MOTIFS
        //
        /////////////////////////////////////////////


        // This has been a bit of a headache: Background_sequences has 4 elements
        // And I need to combine them 1 to 1 with the module and random sequences
        // So each module from_X_to_Y is executed with its correcponding background
        // By using combine with the key the from-to we have one channel with the
        // correct combination. 
        // Peak motifs recives a 5-element tuple in consequence so this can be 
        // used directly as input
        module_with_bg = regulon_sequences.combine(background_sequences, by: 2)
        random_with_bg = random_sequences.combine(background_sequences, by: 2)


        dpm = PEAK_MOTIFS(module_with_bg ,"footDB", "/packages/rsat/public_html/motif_databases/footprintDB/footprintDB.plants.motif.tf")
        rpm = PEAK_MOTIFS_RANDOM(random_with_bg ,"footDB", "/packages/rsat/public_html/motif_databases/footprintDB/footprintDB.plants.motif.tf")

        // We turn the results of PEAK_MOTIFS to MATRIX_SCAN
        // We combine them to be given as input to MATRIX_SCAN, and 
        // transpose them: PEAK_MOTIFS returns all 5 matrices in an array 
        // but each needs their own execution. Transpose turns the channel
        // from [a,b, [c*n times], d] to n elements of shape [a,b,c',d]
        sequences = dpm.combine(background_sequences.seqs, by: 2).transpose()
        r_sequences = rpm.combine(background_sequences.seqs, by: 2).transpose()

        MATRIX_SCAN(sequences, params.pvalue )
        MATRIX_SCAN_RANDOM(r_sequences, params.pvalue)
    }

}