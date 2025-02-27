/*
 * Download MGnify mapseq DB
*/
process GET_MAPSEQ_DB {
    tag "${fn}"
    label 'process_single'
    container 'quay.io/biocontainers/gnu-wget:1.18--hb829ee6_10'

    publishDir "${params.databases.cache_path}", mode: 'copy'
    errorStrategy 'retry'

    input:
    tuple val(meta), val(remote_path), val(md5_path)

    output:
    tuple val(meta), path("${meta.name}"), emit: db_dir

    script:
    fn = remote_path.tokenize('/').last()
    checksum_cmd = ''
    if(!md5_path=='') { 
        checksum_cmd = """
                       wget "${md5_path}" -O md5.2 
                       cat ${fn} | md5sum > md5.1
                       if ! cmp --silent -- md5.1 md5.2; then
                           exit 5
                       fi             
                       # rm md5.*
                       """
    }
    """
    wget "${remote_path}"
    ${checksum_cmd}
    tar -xvzf "${fn}"
    rm "${fn}"
    """
}
