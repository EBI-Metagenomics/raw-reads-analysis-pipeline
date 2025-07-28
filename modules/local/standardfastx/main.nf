
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

    def out_fn = in_f1.getName()
    if (out_fn.endsWith('.gz')) {
        out_fn = out_fn.tokenize('.')[0..<-1].join('.')
    }

    def out_dir = "standardised"

    def reads_cmd = (in_f2 == null) ? "-1 \"${in_f1}\"" : "-1 \"${in_f1}\" -2 \"${in_f2}\""
    def pe_cmd = meta.single_end ? "" : "-p"

    def script = file("${moduleDir}/bin/standardise_fastx.py")

    """
    mkdir -p ${out_dir}
    python ${script} ${args} ${reads_cmd} -o "${out_dir}/${out_fn}" -c ${pe_cmd}
    ls ${out_dir}/* | xargs gzip

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

    def out_fn1 = in_f1.getName()
    def out_fn2 = in_f2.getName()
    def out_dir = "standardised"

    """
    mkdir -p ${out_dir}
    touch "${out_dir}/${out_fn1}"
    touch "${out_dir}/${out_fn2}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        standardfastx: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """
}
