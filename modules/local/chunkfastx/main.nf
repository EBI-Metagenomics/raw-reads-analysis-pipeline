process CHUNKFASTX {
    tag "${meta.id}"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.13'
        : 'quay.io/biocontainers/python:3.13'}"

    input:
    tuple val(meta), path(reads)
    path script

    output:
    tuple val(meta), path("chunked/*"), emit: reads
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def in_f1 = null
    def in_f2 = null
    if (meta.single_end) {
        in_f1 = reads
    } else {
        (in_f1, in_f2) = reads
    }

    def prefix = in_f1.getName().tokenize('.')[0]
    def extension = in_f1.getName().tokenize('.')[1..-1].join('.')
    if (extension.endsWith('.gz')) {
    	extension = extension.tokenize('.')[0..-2].join('.')
    }
    def out_fn = "${prefix}.${extension}"

    def reads_cmd = meta.single_end ? "-1 \"${in_f1}\"" : "-1 \"${in_f1}\" -2 \"${in_f2}\""

    """
    mkdir chunked
    python ${script} ${args} ${reads_cmd} -o "chunked/${out_fn}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chunkfastx: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """

    stub:
    def in_f1 = null
    def in_f2 = null
    if (meta.single_end) {
        in_f1 = reads
    } else {
        (in_f1, in_f2) = reads
    }

    def prefix = in_f1.getName().tokenize('.')[0]
    def extension = in_f1.getName().tokenize('.')[1..-1].join('.')
    def out_fn = "${prefix}.${extension}"
    """
    mkdir chunked
    touch "chunked/${prefix}.${extension}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        renamepairedfastxheaders: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """
}
