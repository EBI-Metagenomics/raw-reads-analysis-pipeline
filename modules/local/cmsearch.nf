/*
 * Infernal cmsearch 1.1.4
*/

process CMSEARCH {

    label 'cmsearch'
    tag "${meta.id}"
    container 'quay.io/biocontainers/infernal:1.1.4--pl5321hec16e2b_1'

    input:
    tuple val(meta), path(sequences)
    tuple val(meta_db), val(covariance_model_database)

    output:
    tuple val(meta), path("${sequences.baseName}*.cmsearch_matches.tbl")

    script:
    """
    cmsearch \
    --cpu ${task.cpus} \
    --cut_ga \
    --noali \
    --hmmonly \
    -Z 1000 \
    -o /dev/null \
    --tblout ${sequences.baseName}.cmsearch_matches.tbl \
    ${covariance_model_database} \
    ${sequences}
    """
}
