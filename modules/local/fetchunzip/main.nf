process FETCHUNZIP {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--hb829ee6_10'
        : 'quay.io/biocontainers/gnu-wget:1.18--hb829ee6_10'}"

    publishDir "${params.databases.cache_path}", mode: 'copy'
    errorStrategy 'retry'

    input:
    tuple val(meta), val(db_map), path(fp)

    output:
    tuple val(meta), path("${meta.id}")

    script:
    def fn = fp.tokenize('/').last()
    if (fn[-7..-1]=='.tar.gz') {
        """
        #!/bin/bash
        # wget "${db_map.remote_path}"
        tar -xvzf "${fp}"
        # rm "${fn}"
        """
    } else {
        """
        #!/bin/bash
        cp "${fp}" .
        """
    }
}
