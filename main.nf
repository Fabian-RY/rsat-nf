#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC } from './modules/nf-core/fastqc/main.nf'
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main.nf'


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

    // Workflow steps:
    // FASTQC: Performs basic QC of the RNA samples. This step is totally independent, and can be executed inmediatelly at beggining
    // Then Trimmomatic starts QC of the samples
    preprocess = FASTQC(samplesheet_ch)
    Trimmomatic_result = TRIMMOMATIC(samplesheet_ch)
    postprocess = FASTQ(Trimmomatic_result.out.trimmed_reads)

}