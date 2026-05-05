#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC } from './modules/nf-core/fastqc/main.nf'
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main.nf'


// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

    // input channels
    samplesheet_ch = Channel
        .fromPath(params.samplesheet)
        .ifEmpty { error "No samples found in: ${params.samplesheet}" }
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample_id.trim(), single_end: true]
            def fastq_file = file(row.fastq.trim())
            [meta, fastq_file]
        }

    // Split channel for parallel execution
//    samplesheet_ch.into { ch_fastqc; ch_trimmomatic }

    // Workflow steps
    FASTQC(samplesheet_ch)
    TRIMMOMATIC(samplesheet_ch)

    publish:
    output = FASTQC.out.html.mix(FASTQC.out.zip).mix(TRIMMOMATIC.out.trimmed_reads).mix(TRIMMOMATIC.out.summary)
}
