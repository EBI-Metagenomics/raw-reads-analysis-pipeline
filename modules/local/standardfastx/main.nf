
process STANDARDFASTX {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'
        : 'biocontainers/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("standardised/*"), emit: reads
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def in_f1 = null
    def in_f2 = null
    if (reads.size() == 1) {
        in_f1 = reads[0]
    }
    else {
        (in_f1, in_f2) = reads
    }

    def prefix = in_f1.getName().tokenize('.')[0]
    def extension = in_f1.getName().tokenize('.')[1..-1].join('.')
    if (extension.endsWith('.gz')) {
        extension = extension.tokenize('.')[0..-2].join('.')
    }
    def out_fn = "${prefix}.${extension}"

    def reads_cmd = meta.single_end ? "-1 \"${in_f1}\"" : "-1 \"${in_f1}\" -2 \"${in_f2}\""
    def pe_cmd = meta.single_end ? "" : "-p"

    def script = file("${moduleDir}/bin/standardise_fastx.py")

    """
    mkdir standardised
    python ${script} ${args} ${reads_cmd} -o "standardised/${out_fn} -c ${pe_cmd}"
    ls standardised | xargs gzip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        standardfastx: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """

    stub:
    def in_f1 = null
    def in_f2 = null
    if (reads.size() == 1) {
        in_f1 = reads
    }
    else {
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
        standardfastx: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """
}
