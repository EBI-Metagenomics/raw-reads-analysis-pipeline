process GZIPALL {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(files, name: 'in_files/*')

    output:
    tuple val(meta), path('out_files/*'), emit: files
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    
    """
    mkdir out_files
    for f in in_files/* ; do gzip -c \"\$(readlink "\$f")\" > "out_files/\$(basename \"\$f\").gz" ; done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzipall: \$(gzip --version |& sed '1!d ; s/gzip //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    mkdir out_files
    for f in in_files/* ; do touch "out_files/\$f.gz" ; done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzipall: \$(gzip --version |& sed '1!d ; s/gzip //')
    END_VERSIONS
    """
}
