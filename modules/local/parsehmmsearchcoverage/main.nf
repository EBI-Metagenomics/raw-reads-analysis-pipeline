process PARSEHMMSEARCHCOVERAGE {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'
        : 'biocontainers/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(domtbl_file)

    output:
    tuple val(meta), path(out_fp), emit: tsv
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    out_fp = "${meta.id}_pfam_coverage.tsv"
    """
    gunzip -c ${domtbl_file} | python hmmer_domtbl_parse_coverage.py ${args} -i - -o ${out_fp}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsehmmsearchcoverage: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    out_fp = "${meta.id}_pfam_coverage.tsv"

    """
    touch ${out_fp}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsehmmsearchcoverage: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
