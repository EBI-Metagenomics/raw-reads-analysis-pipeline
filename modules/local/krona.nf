/*
 * Krona 2.7.1
*/

process KRONA {

    label 'krona'
    tag "${meta.seq_id} - ${meta.db_id}"
    
    publishDir(
        "${params.outdir}/taxonomy/${meta.seq_id}/${meta.db_id}",
        mode: 'copy',
        pattern: "*krona.html"
    )

    container 'quay.io/biocontainers/krona:2.7.1--pl5321hdfd78af_7'

    input:
    tuple val(meta), path(tsv), path(txt), path(tsv_notaxid)

    output:
    tuple val(meta), path("*krona.html")

    script:
    """
    ktImportText -o "${meta.seq_id}_${meta.db_id}_krona.html" $tsv
    """
}
