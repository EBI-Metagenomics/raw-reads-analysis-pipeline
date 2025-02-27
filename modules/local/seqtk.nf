/*
 * seqtk 1.3
*/

process SEQTK {
    publishDir "${params.outdir}/", mode: 'copy'
    label 'seqtk'
    tag "${meta.id}"
    container 'quay.io/biocontainers/seqtk:1.3--h7132678_4'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("${meta.id}.fasta")

    script:
    """
    seqtk seq -a ${reads} > ${meta.id}.fasta
    """
}
