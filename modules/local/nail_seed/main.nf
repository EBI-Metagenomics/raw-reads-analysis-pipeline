process NAIL_SEED {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/nail:0.7.1--hab7d0fd_0'
        : 'biocontainers/nail:0.7.1--hab7d0fd_0'}"

    input:
    tuple val(meta), path(orfs), path(hmm)

    output:
    tuple val(meta), path("${prefix}.seeds.tsv"), emit: seeds
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def target = orfs.name.endsWith('.gz') ? orfs.name[0..-4] : orfs.name
    def unzip = orfs.name.endsWith('.gz') ? "gunzip -c ${orfs} > ${target}" : ''
    """
    ${unzip}

    nail search \\
        --only-seed \\
        ${args} \\
        -t ${task.cpus} \\
        --tmp-dir tmp \\
        --seeds-out ${prefix}.seeds.tsv \\
        ${hmm} \\
        ${target}

    rm -rf tmp ${orfs.name.endsWith('.gz') ? target : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nail: \$(nail --version | sed 's/nail //')
        mmseqs: \$(mmseqs version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.seeds.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nail: stub
        mmseqs: stub
    END_VERSIONS
    """
}
