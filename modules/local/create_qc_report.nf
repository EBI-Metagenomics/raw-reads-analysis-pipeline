/*
 * fastp json parser to report number of reads after every filtering
*/
process QC_REPORT {

    publishDir "${params.outdir}/qc", mode: 'copy'

    container 'quay.io/biocontainers/python:3.11'
    tag "${meta.id}"
    label 'report'

    input:
    tuple val(meta), val(reports)

    output:
    tuple val(meta), path("qc_summary")

    script:
    def inputs = ""
    if (reports.samples.mode == "paired") {
        inputs = " --overlap-counts ${reports.merged}" 
    }
    else {
        inputs = ""
    }

    """
    collect_counts.py \
    --qc-json ${reports.fastp} \
    --decontamination-counts ${reports.decontam} \
    ${inputs} -o "qc_summary"
    """
}
