#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateParameters; paramsHelp; paramsSummaryLog; } from 'plugin/nf-schema'
include { PIPELINE } from './workflows/rrap.nf'

workflow {
    // Validate input parameters
    validateParameters()

    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    PIPELINE ()
}
