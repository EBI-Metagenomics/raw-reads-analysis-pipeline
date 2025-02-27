/*
 * MAPseq 2.1.1
*/
process MAPSEQ {

    publishDir(
        "${params.outdir}/taxonomy/${meta.seq_id}/${meta.db_id}",
        mode: 'copy',
        pattern: "*.mseq*"
    )

    container 'quay.io/biocontainers/mapseq:2.1.1--ha34dc8c_0'

    label 'mapseq'
    tag "${meta.seq_id}"
    
    input:
    tuple val(meta),val(sequence), path(db_dir)

    output:
    tuple val(meta), path("${sequence.baseName}.mseq")

    script:
    fasta_db = "${db_dir}/${params.databases[meta.db_id].files.fasta}"
    tax_db = "${db_dir}/${params.databases[meta.db_id].files.tax}"
    """
    mapseq \
        ${sequence} \
        ${fasta_db} \
        ${tax_db} \
        -nthreads ${task.cpus} \
        -tophits 80 \
        -topotus 40 \
        -outfmt 'simple' > ${sequence.baseName}.mseq
    """
}


