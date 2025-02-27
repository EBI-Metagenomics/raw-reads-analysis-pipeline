process FETCH_URL {
    tag "$meta.id"

    label 'process_single'
    container 'quay.io/biocontainers/gnu-wget:1.18--hb829ee6_10'

    publishDir "${params.reads_cache_path}", mode: 'copy'
    errorStrategy 'retry'

    input:
    tuple val(meta), val(url), val(md5), val(read_type)

    output:
    tuple val(meta), path("${meta.id}/${fn}"), val(read_type)

    script:
    fn = url.tokenize('/').last()
    checksum_cmd = ''
    // if(!(md5=='')) {
    //     checksum_cmd = """
    //                    if \$dl_md5 -nq $md5; then
    //                        exit 5
    //                    fi
    //                    """
    // }
    """
    mkdir -p ${meta.id}
    wget ${url} -O ${meta.id}/${fn}
    # dl_md5="\$(cat $meta.id/$fn | md5sum | awk '{print \$1}')"
    ${checksum_cmd}
    """
}
