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
include { RETRIEVE_ALL_SEQUENCES } from './modules/local/rsat-retrieve-sequences/main.nf'
include { PURGE_SEQUENCES } from './modules/local/rsat-purge-sequences/main.nf'
include { DYAD } from './modules/local/rsat-dyad-analysis/main.nf'
include { OLIGO } from './modules/local/rsat-oligo-analysis/main.nf'
include { PEAK_MOTIFS as PEAK_MOTIFS_DYAD} from './modules/local/rsat-peak-motifs/main.nf'
include { PEAK_MOTIFS as PEAK_MOTIFS_OLIGO} from './modules/local/rsat-peak-motifs/main.nf'
include { MATRIX_SCAN as MATRIX_SCAN_DYAD } from './modules/local/rsat-matrix-scan/main.nf'
include { MATRIX_SCAN as MATRIX_SCAN_OLIGO } from './modules/local/rsat-matrix-scan/main.nf'



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
        // Random Gene selection without download
        n_clusters_ch = Channel
                        .from(1..params.clusters)
                        .each {
                            i -> "random_cluster_${1}"
                        }
        all_sequences = RETRIEVE_ALL_SEQUENCES(params.organism, params.retrieve_seq_output,
                params.feature_type, params.retrieve_seq_type, params.retrieve_seq_format, params.retrieve_seq_label, params.retrieve_seq_from,
                params.retrieve_seq_to)
        clusters = RANDOM_GENES(n_clusters_ch, params.number, params.organism, params.feature_type)         
        module_sequences = RETRIEVE_RANDOM_SEQUENCES(params.organism, clusters.outfile, params.retrieve_seq_output,
                params.feature_type, params.retrieve_seq_type, params.retrieve_seq_format, params.retrieve_seq_label, params.retrieve_seq_from,
                params.retrieve_seq_to)
        sequences = RETRIEVE_SEQUENCES(params.organism, modules_ch, params.retrieve_seq_output,
                params.feature_type, params.retrieve_seq_type, params.retrieve_seq_format, params.retrieve_seq_label, params.retrieve_seq_from,
                params.retrieve_seq_to)
        purged_sequences = PURGE_SEQUENCES(sequences, params.retrieve_seq_format)
        dyad_ch = DYAD(purged_sequences, params.organism)
        oligo_ch = OLIGO(purged_sequences, params.organism)
        pk_dyad_ch = PEAK_MOTIFS_DYAD(dyad_ch, module_sequences)
        pk_oligo_ch = PEAK_MOTIFS_OLIGO(oligo_ch, module_sequences)
        MATRIX_SCAN_DYAD(pk_dyad_ch, params.organism)
        MATRIX_SCAN_OLIGO(pk_oligo_ch, params.organism)
    }


    if(params.preprocess_samples){

        // input channels
        samplesheet_ch = Channel.fromPath(params.samplesheet)
            .ifEmpty { error "Error: No samples found in: ${params.samplesheet}"}
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.sample_id.trim(), single_end: true]
                def fastq_file = file(row.fastq.trim())
                [meta, fastq_file]
            }

        kallisto_index_ch = Channel.of(tuple(params.organism, params.transcriptome_fa))
        kallisto_quant_gtf = Channel.of(params.transcriptome_gtf)
        kallisto_quant_insert_length_ch = Channel.of(params.fragment_length)
        kallisto_quant_insert_sd_ch = Channel.of(params.fragment_sd)

        // Workflow steps:
        // Input: Fasta files of samples
        // FASTQC: Performs basic QC of the RNA samples. This step is totally independent, and can be executed inmediatelly at beggining
        // Then Trimmomatic starts QC of the samples
        // Then FastQC of the remaining reads
        preprocess = FASTQC_RAW_READS(samplesheet_ch)
        Trimmomatic_result = TRIMMOMATIC(samplesheet_ch)
        postprocess = FASTQC_TRIMMERED_READS(Trimmomatic_result.trimmed_reads)

        // Input: Transcriptome indexing
        // Transcriptome must be indexed before running the quantification
        k_index = KALLISTO_INDEX(kallisto_index_ch)
        
        // Quantification: With reads and index we can start quantification
        // TO DO: Implement chromosomes as a list of inputs from a file or sth
        KALLISTO_QUANT(Trimmomatic_result.trimmed_reads, k_index.index, kallisto_quant_gtf , [], kallisto_quant_insert_length_ch, kallisto_quant_insert_sd_ch )

    }
    
    if(params.preprocess_samples){

        // MULTIQC Final Report
        ch_multiqc_files = Channel.empty()
            .mix(
                FASTQC_RAW_READS.out.zip.map { meta, zip -> zip },
                FASTQC_RAW_READS.out.html.map { meta, html -> html },
                FASTQC_TRIMMERED_READS.out.zip.map { meta, zip -> zip },
                FASTQC_TRIMMERED_READS.out.html.map { meta, html -> html },
                TRIMMOMATIC.out.trim_log.map { meta, log -> log },
                TRIMMOMATIC.out.out_log.map { meta, log -> log },
                TRIMMOMATIC.out.summary.map { meta, summary -> summary },
                KALLISTO_QUANT.out.log.map { meta, log -> log },
                KALLISTO_QUANT.out.json_info.map { meta, json -> json }
            )
            .collect()
            .map { files -> [[id:'multiqc'], files, [], [], [], []] }
        MULTIQC(ch_multiqc_files)
    }
}