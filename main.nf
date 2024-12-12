#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RNASEQ } from './workflows/rnaseq'

workflow {
    RNASEQ()
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Execution duration: $workflow.duration"
}