process MULTIQC {

    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'
    tag "${meta.id}"
    container 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'

    input:
    tuple val(meta), path(fastp_json), path(motus_log)

    output:
    tuple val(meta), path("multiqc_report.html"), path("multiqc_data")

    script:
    """
    multiqc --module fastp --module motus .
    """
}
