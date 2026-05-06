#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW_READS} from './modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_TRIMMERED_READS  } from './modules/nf-core/fastqc/main.nf'
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main.nf'
include { KALLISTO_INDEX } from './modules/nf-core/kallisto/index/main.nf' 
include { KALLISTO_QUANT } from './modules/nf-core/kallisto/quant/main.nf'

// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

    // input channels
    samplesheet_ch = Channel.fromPath(params.samplesheet)
        .ifEmpty { error "Error: No samples found in: ${params.samplesheet}"}
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample_id.trim(), single_end: true]
            def fastq_file = file(row.fastq.trim())
            [meta, fastq_file]
        }

    kallisto_index_ch = Channel.of(tuple(params.transcriptome_name, params.transcriptome))

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

}