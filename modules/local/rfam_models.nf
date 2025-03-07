/*
 * Download MGnify Rfam DB
*/
process GET_RFAM_DB {
    
    tag "${meta.id}"
    label 'get_rfam_db'

    container 'quay.io/biocontainers/gnu-wget:1.18--hb829ee6_10'
    publishDir "${params.databases.cache_path}", mode: 'copy'

    input:
    tuple val(meta), val(remote_path), val(md5_path)

    output:
    tuple val(meta), path("${params.databases.rfam.name}"), emit: db_dir

    script:
    fn = remote_path.tokenize('/').last()
    """
    wget "${remote_path}" 
    tar -xvzf "${fn}"
    rm "${fn}"
    """
}
