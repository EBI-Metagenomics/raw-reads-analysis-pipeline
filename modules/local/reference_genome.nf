/*
 * Download reference genome, default: HG38
*/
process GET_REFERENCE_GENOME {

    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/gnu-wget:1.18--hb829ee6_10'

    publishDir "${params.databases.cache_path}", mode: 'copy'
    errorStrategy 'retry'

    input:
    tuple val(meta), val(remote_path), val(md5_path)

    output:
    tuple val(meta), path("${params.databases.host_genome.name}")

    script:
    def fn = remote_path.tokenize('/').last()
    def checksum_cmd = ''
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
