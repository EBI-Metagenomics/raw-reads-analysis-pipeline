process ADDHEADER {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0' :
        'biocontainers/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(in_fp)
    val header_str

    output:
    tuple val(meta), path("output/*"), emit: file_with_header

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    fn=\$(basename $in_fp)
    mkdir output
    out_fp=output/\$fn
    echo -e "$header_str" > \$out_fp
    cat $in_fp >> \$out_fp
    """

    stub:
    def args = task.ext.args ?: ''

    """
    fn=\$(basename $in_fp)
    mkdir output
    out_fp = output/\$fn
    touch \$out_fp
    """
}
