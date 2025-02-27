/*
 * mapseq2biom python script converter v1.0.0
*/

process MAPSEQ2BIOM {

    publishDir(
        path: "${params.outdir}/taxonomy/${meta.seq_id}/${otu_label}",
        pattern: "${mapseq_out.baseName}.*",
        mode: 'copy'
    )

    container 'quay.io/biocontainers/python:3.11'

    label 'mapseq2biom'
    tag "${meta.seq_id}"
    
    input:
        tuple val(meta), path(mapseq_out), val(sequence), val(db_dir)
    
    output:
        tuple val(meta), path("${mapseq_out.baseName}.tsv"), path("${mapseq_out.baseName}.txt"), path("${mapseq_out.baseName}.notaxid.tsv")

    script:
    otu_label = meta.db_id
    otu_ref = "${db_dir}/${params.databases[meta.db_id].files.otu}"
    """
    mapseq2biom.py \
        --out-file ${mapseq_out.baseName}.tsv \
        --krona ${mapseq_out.baseName}.txt \
        --no-tax-id-file ${mapseq_out.baseName}.notaxid.tsv \
        --taxid \
        --label ${otu_label} \
        --query ${mapseq_out} \
        --otu-table ${otu_ref}
    """
}
