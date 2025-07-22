process FASTAEMBEDLENGTH {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'
        : 'biocontainers/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${out_fn}.gz"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def base = fasta.getName().tokenize('.')[0..-2].join('.')
    // def ext = fasta.getName().tokenize('.')[-1]
    out_fn = "${base}.renamed"

    def script = file("${moduleDir}/bin/fastx_embed_length.py")

    """
    python ${script} -i ${fasta} -o ${out_fn} --no_output_gzip ${args}
    gzip ${out_fn}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """

    stub:
    def base = fasta.getName().tokenize('.')[0..-2].join('.')
    out_fn = "${base}.renamed"

    """
    touch ${out_fn}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """
}
