/*
 * Extract models from cmsearched fasta
*/

process EXTRACT_MODELS {

    publishDir "${params.outdir}/cmsearch/", mode:'copy'
    
    label 'extract_coords'
    tag "${meta.id}"
    container 'quay.io/biocontainers/biopython:1.75'

    input:
    tuple val(meta), path(sequences) 

    output:
    tuple val(meta), path("sequence-categorisation/"), path("sequence-categorisation/${meta.id}_SSU.fasta"), path("sequence-categorisation/${meta.id}_LSU.fasta")

    script:
    """
    get_subunits.py -i ${sequences} -n $meta.id
    """
}
