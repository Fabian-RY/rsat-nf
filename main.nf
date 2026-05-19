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


// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {


    // Checks to download the organism or not
    if (params.download_organism){
        download_result = DOWNLOAD_ORGANISM(params.server, params.organism)
        def organism_downloaded = download_result.map {org, server -> org}
    }
    // Random Gene selection without download
    genes = RANDOM_GENES(params.number, params.organism, params.feature_type)
    

    sequences = RETRIEVE_SEQUENCES

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