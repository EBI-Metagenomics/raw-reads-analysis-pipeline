process FETCHUNZIP {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--hb829ee6_10':
        'biocontainers/gnu-wget:1.18--hb829ee6_10' }"

    publishDir "${params.databases.cache_path}", mode: 'copy'
    errorStrategy 'retry'
    
    input:
    tuple val(meta), val(db_map)
    
    output:
    tuple val(meta), path("${db_map.name}") 

    script:
    def fn = db_map.remote_path.tokenize('/').last()
    """
    #!/bin/bash
    wget "${db_map.remote_path}"
    tar -xvzf "${fn}"
    rm "${fn}"
    """
}
