/*
 * qc_stats for fasta or fastq file (was named MGRAST_base.py)
*/
process QC_STATS {

    publishDir "${params.outdir}/qc", mode: 'copy'

    container 'quay.io/biocontainers/biopython:1.75'

    label 'qc_summary'
    tag "${meta.id}"

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("statistics")

    script:
    """
    qc_summary.py -i ${sequence} --output-dir "statistics"
    """
}
