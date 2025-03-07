#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PIPELINE } from './workflows/rrap.nf'

workflow {
    PIPELINE ()
}
