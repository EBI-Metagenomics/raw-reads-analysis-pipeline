process ADDHEADER_GZIP {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0' :
        'biocontainers/mgnify-pipelines-toolkit:0.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(in_fp), val(out_fn)
    val header_str
    val gzip

    output:
    tuple val(meta), path("output/*"), emit: file_with_header

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    if (!out_fn) {
        if (in_fp) {
            out_fn = in_fp.name
        }
        else {
            out_fn = "file"
        }
    }

    def echo_cmd = "touch \$out_fp"
    if ((header_str instanceof String) && (header_str.length()>0)){
        echo_cmd = "echo -e \"${header_str}\" > \$out_fp"
    }

    def cat_cmd = ""
    if (in_fp) {
        cat_cmd = "cat ${in_fp} >> \$out_fp"
    }

    def gzip_cmd = ""
    if (gzip) {
        gzip_cmd = "gzip \$out_fp"
    } 

    """
    mkdir -p output
    out_fp=output/$out_fn
    ${echo_cmd}
    ${cat_cmd}
    ${gzip_cmd}
    """

    stub:
    def args = task.ext.args ?: ''

    if (!out_fn) {
        out_fn = in_fp.name
    }
    if (gzip) {
        out_fn = out_fn + '.gz'
    }

    """
    mkdir -p output
    out_fp = output/$out_fn
    touch \$out_fp
    """
}
