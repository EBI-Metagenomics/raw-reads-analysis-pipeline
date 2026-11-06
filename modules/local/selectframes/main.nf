process SELECTFRAMES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'
        : 'biocontainers/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(reads), path(orfs), path(calls)

    output:
    tuple val(meta), path("${prefix}_selected.faa.gz"), emit: orfs
    tuple val(meta), path("${prefix}_selectframes.json"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def script = file("${moduleDir}/bin/select_frames.py")
    """
    python ${script} \\
        ${args} \\
        --reads ${reads} \\
        --orfs ${orfs} \\
        --calls ${calls} \\
        --out ${prefix}_selected.faa.gz \\
        --stats ${prefix}_selectframes.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -n | gzip > ${prefix}_selected.faa.gz
    echo '{}' > ${prefix}_selectframes.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
