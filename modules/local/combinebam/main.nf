process COMBINEBAM {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0'
        : 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'}"

    input:
    tuple val(meta), path(input_files, stageAs: "inputs/?/*")

    output:
    tuple val(meta), path("${output_name}"), emit: concatenated_result
    path "versions.yml"                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    output_name = (input_files instanceof List ? input_files[0] : input_files).name.tokenize('/')[-1]

    """
    samtools cat -o ${output_name} \$(ls inputs/*/*) ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combinebam: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    output_name = input_files[0].name

    """
    touch ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combinebam: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
