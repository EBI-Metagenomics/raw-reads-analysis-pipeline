process MOTUS2KRONA {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'
        : 'biocontainers/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(motus_out)

    output:
    tuple val(meta), path("*_krona.tsv"), emit: krona
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def script = file("${moduleDir}/bin/motus2krona.py")

    """
    python ${script} \\
        ${args} \\
        -i "${motus_out}" \\
        -o "${meta.id}_krona.tsv" \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version |& sed '1!d ; s/Python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch "${meta.id}_krona.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version |& sed '1!d ; s/Python //')
    END_VERSIONS
    """
}
