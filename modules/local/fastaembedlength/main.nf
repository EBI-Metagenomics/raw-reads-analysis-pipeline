process FASTAEMBEDLENGTH {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${out_fn}"), emit: fasta
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def base = fasta.getName().tokenize('.')[0..-2].join('.')
    def ext = fasta.getName().tokenize('.')[-1]
    out_fn = "${base}.renamed.${ext}"

    """
    awk '/^>/ { if (name) {printf("%s_length=%d\n%s", name, len, seq)} name=$0; seq=""; len = 0; next}
    NF > 0 {seq = seq $0 "\n"; len += length()}
    END { if (name) {printf("%s_length=%d\n%s", name, len, seq)} }' ${fasta} > ${out_fn}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args = task.ext.args ?: ''
    def base = fasta.getName().tokenize('.')[0..-2].join('.')
    def ext = fasta.getName().tokenize('.')[-1]
    out_fn = "${base}.renamed.${ext}"
    
    """
    touch ${out_fn}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """
}
